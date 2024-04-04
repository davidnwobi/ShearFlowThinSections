import pandas as pd
import pytest
from shape import Shape
from sympy import symbols, expand, integrate
from data_setup import dimensions
from open_section_solver import OpenSectionSolver
from closed_section_solver import SingleClosedSectionSolver


def test_ex1_closed_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 2 * a, 'h': 2 * a})
    shape = Shape()
    solver = SingleClosedSectionSolver(shape)
    solver.solve()

    S01 = shape.elements[(0, 1)].S
    S12 = shape.elements[(1, 2)].S
    S23 = shape.elements[(2, 3)].S
    S34 = shape.elements[(3, 4)].S
    S45 = shape.elements[(4, 5)].S
    S_z = dimensions['S_z']

    assert expand(shape.elements[(0, 1)].qb.subs(S01, 0) - 0) == 0
    assert expand(shape.elements[(0, 1)].qb.subs(S01, a) + 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(1, 2)].qb.subs(S12, 0) + 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(1, 2)].qb.subs(S12, 2 * a) + 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(2, 3)].qb.subs(S23, 0) + 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(2, 3)].qb.subs(S23, 2 * a) - 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(3, 4)].qb.subs(S34, 0) - 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(3, 4)].qb.subs(S34, 2 * a) - 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(4, 5)].qb.subs(S45, 0) - 3 * S_z / (16 * a)) == 0
    assert expand(shape.elements[(4, 5)].qb.subs(S45, a) - 0) == 0

    assert (shape.qo - S_z / (16 * a)) == 0
