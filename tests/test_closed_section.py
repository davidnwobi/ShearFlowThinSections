from shape import Shape
from sympy import symbols, expand, sqrt, nsimplify
from data_setup import dimensions
from closed_section_solver import SingleClosedSectionSolver


def test_ex1_closed_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 2 * a, 'h': 2 * a})
    shape = Shape(nodes='ex1_closed_section/Nodes.xlsx',
                  elements='ex1_closed_section/Elements.xlsx',
                  forces='ex1_closed_section/Forces.xlsx',
                  sections='ex1_closed_section/ClosedSections.xlsx',
                  shear_center_ref='ex1_closed_section/Shear Center Ref.xlsx')

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

def test_ex2_closed_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 2 * a, 'h': sqrt(3) * a})
    shape = Shape(nodes='ex2_closed_section/Nodes.xlsx',
                  elements='ex2_closed_section/Elements.xlsx',
                  forces='ex2_closed_section/Forces.xlsx',
                  sections='ex2_closed_section/ClosedSections.xlsx',
                  shear_center_ref='ex2_closed_section/Shear Center Ref.xlsx')

    solver = SingleClosedSectionSolver(shape)
    solver.solve()

    S_y = dimensions['S_y']
    assert expand(shape.elements[(1, 2)].Qb - S_y / 18) == 0
    assert expand(shape.elements[(0, 2)].Qb - 5 * S_y / 18) == 0
    assert expand(shape.elements[(2, 3)].Qb - 7 / 9 * S_y) == 0
    assert expand(shape.elements[(3, 5)].Qb - 5 / 18 * S_y) == 0
    assert expand(shape.elements[(3, 4)].Qb - S_y / 18) == 0

    assert expand(nsimplify(shape.qo, rational=True) - 7 / 9 * (S_y / a)) == 0
