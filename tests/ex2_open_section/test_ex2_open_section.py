import pytest
from shape import Shape
from sympy import symbols, expand, integrate
from data_setup import dimensions
from open_section_solver import OpenSectionSolver


def test_ex2_open_section():
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

    assert expand(shape.elements[(0, 1)].tz - (-2 * a * t)) == 0
    assert expand(shape.elements[(0, 1)].int_tz - (-2 * a * t * S01)) == 0
    assert shape.elements[(0,1)].qs.subs(S01, a) == 3*Sz/(7*a)
    assert shape.elements[(0,1)].Q == 3*Sz/14

    S13 = shape.elements[(1, 3)].S
    l13 = shape.elements[(1, 3)].length

    assert expand(shape.elements[(1, 3)].tz - t * (S13 - a)) == 0
    assert expand(shape.elements[(1, 3)].int_tz - (t * (S13 ** 2 / 2 - a * S13) - 2 * a ** 2 * t)) == 0
    assert expand(shape.elements[(1, 3)].qs.subs(S13, a) - 15 * Sz / (28 * a)) == 0
    assert expand(shape.elements[(1, 3)].Q - Sz/2) == 0

    S23 = shape.elements[(2, 3)].S
    l23 = shape.elements[(2, 3)].length

    assert expand(shape.elements[(2, 3)].tz - 0) == 0
    assert expand(shape.elements[(2, 3)].int_tz - 0) == 0

    S34 = shape.elements[(3, 4)].S
    l34 = shape.elements[(3, 4)].length

    assert expand(shape.elements[(3, 4)].tz - S34 * t) == 0
    assert expand(shape.elements[(3, 4)].int_tz - (S34 ** 2 / 2 * t - 5 * a ** 2 * t / 2)) == 0
    assert expand(shape.elements[(3, 4)].qs.subs(S34, 0) - 15 * Sz / (28 * a)) == 0
    assert expand(shape.elements[(3, 4)].Q - Sz/2) == 0

    S45 = shape.elements[(4, 5)].S
    l45 = shape.elements[(4, 5)].length

    assert expand(shape.elements[(4, 5)].tz - 2 * a * t) == 0
    assert expand(shape.elements[(4, 5)].int_tz - 2 * a * t * (S45 - a)) == 0
    assert expand(shape.elements[(4, 5)].qs.subs(S45, 0) - 3*Sz/(7*a)) == 0
    assert expand(shape.elements[(4, 5)].Q - 3*Sz/14) == 0


