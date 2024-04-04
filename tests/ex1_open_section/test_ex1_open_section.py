import pytest
from shape import Shape
from sympy import symbols, expand, integrate
from data_setup import dimensions
from open_section_solver import OpenSectionSolver


def test_ex1_open_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': a, 'h': 2 * a})
    shape = Shape(nodes='Nodes.xlsx', elements='Elements.xlsx', forces="Forces.xlsx")
    solver = OpenSectionSolver(shape)

    solver.solve()
    t = dimensions['t']
    Sz = dimensions['S_z']
    S01 = shape.elements[(0, 1)].S
    l01 = shape.elements[(0, 1)].length
    assert shape.Iy == 14 * a ** 3 * t / 3
    assert shape.elements[(0,1)].tz == -2*a*t
    assert shape.elements[(0,1)].int_tz == -2*a*t*S01
    assert shape.elements[(0,1)].qs.subs(S01, a) == 3*Sz/(7*a)
    assert shape.elements[(0,1)].Q == 3*Sz/14

    l12 = shape.elements[(1,2)].length
    S12 = shape.elements[(1, 2)].S
    assert expand(shape.elements[(1, 2)].tz - t*(S12-a)) == 0
    assert expand(shape.elements[(1, 2)].int_tz - (t*(S12**2/2-a*S12)-2*a**2*t)) == 0
    assert expand(shape.elements[(1, 2)].qs.subs(S12, a)-15*Sz/(28*a)) == 0
    assert expand(shape.elements[(1, 2)].Q - Sz) == 0

    l23 = shape.elements[(2, 3)].length
    S23 = shape.elements[(2, 3)].S
    assert shape.elements[(2, 3)].tz == 2*a*t
    assert expand(shape.elements[(2, 3)].int_tz - 2*a*t*(S23-a)) == 0
    assert expand(shape.elements[(2, 3)].qs.subs(S23, 0) - 3*Sz/(7*a)) == 0
    assert expand(shape.elements[(2, 3)].Q - 3*Sz/14) == 0


