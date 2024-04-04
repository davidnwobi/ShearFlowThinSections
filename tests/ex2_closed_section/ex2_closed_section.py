import pytest
from shape import Shape
from sympy import symbols, expand, integrate, sqrt, nsimplify
from data_setup import dimensions
from closed_section_solver import SingleClosedSectionSolver


def test_ex2_closed_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 2 * a, 'h': sqrt(3) * a})
    shape = Shape()
    solver = SingleClosedSectionSolver(shape)
    solver.solve()

    S_y = dimensions['S_y']
    assert expand(shape.elements[(1, 2)].Qb - S_y / 18) == 0
    assert expand(shape.elements[(0, 2)].Qb - 5 * S_y / 18) == 0
    assert expand(shape.elements[(2, 3)].Qb - 7 / 9 * S_y) == 0
    assert expand(shape.elements[(3, 5)].Qb - 5 / 18 * S_y) == 0
    assert expand(shape.elements[(3, 4)].Qb - S_y / 18) == 0

    assert expand(nsimplify(shape.qo, rational=True) - 7 / 9 * (S_y / a)) == 0
