import pytest
from shape import Shape
from sympy import symbols, expand, integrate
from data_setup import dimensions
from open_section_solver import OpenSectionSolver


def test_ex3_open_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 2 * a, 'h': 2 * a})
    shape = Shape(nodes='Nodes.xlsx', elements='Elements.xlsx', forces="Forces.xlsx")
    solver = OpenSectionSolver(shape)

    solver.solve()
    t = dimensions['t']
    Sz = dimensions['S_z']
    S01 = shape.elements[(0, 1)].S
    l01 = shape.elements[(0, 1)].length
    assert shape.Iy == 37 * a ** 3 * t / 6
    assert shape.Iz == 15 * a ** 3 * t / 8
    assert shape.Iyz == -11 * a ** 3 * t / 4

    assert expand(shape.elements[(0, 1)].int_ty - 4 * t * (7 / 8 * a * S01 - S01 ** 2 / 2)) == 0
    assert expand(shape.elements[(0, 1)].int_tz - (- 3 * a * t * S01)) == 0
    assert expand(shape.elements[(0, 1)].qs.subs(S01,l01) - 3 * Sz / (8 * a)) == 0

    l12 = shape.elements[(1, 2)].length
    S12 = shape.elements[(1, 2)].S

    assert expand(shape.elements[(1, 2)].int_ty - (-a / 8 * t * S12 + 3 / 2 * a ** 2 * t)) == 0
    assert expand(shape.elements[(1, 2)].int_tz - (t * (-3 * a * S12 / 4 + S12 ** 2 / 2) - 3 * a ** 2 * t)) == 0
    assert expand(shape.elements[(1, 2)].qs.subs(S12, 0) - 3 * Sz / (8 * a)) == 0
    assert expand(shape.elements[(1, 2)].qs.subs(S12, l12) - 5 * Sz / (16 * a)) == 0

    l23 = shape.elements[(2, 3)].length
    S23 = shape.elements[(2, 3)].S

    assert expand(shape.elements[(2, 3)].int_ty - (2 * t * (-a / 8 * S23 - S23 ** 2 / 2) + 5 / 4 * a ** 2 * t)) == 0
    assert expand(shape.elements[(2, 3)].int_tz - (5 / 2 * a * t * S23 - 5 / 2 * a ** 2 * t)) == 0
    assert expand(shape.elements[(2, 3)].qs.subs(S23, 0) - 5 * Sz / (16 * a)) == 0
