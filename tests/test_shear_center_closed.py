from shape import Shape
from sympy import symbols, expand, nsimplify
from data_setup import dimensions
from closed_section_solver import SingleClosedSectionSolver


def test_ex1_shear_center_closed():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 5 * a, 'h': 4.8 * a})
    shape = Shape(
        nodes="ex1_shear_center_closed/Nodes.xlsx",
        elements="ex1_shear_center_closed/Elements.xlsx",
        forces="ex1_shear_center_closed/Forces.xlsx",
        sections="ex1_shear_center_closed/ClosedSections.xlsx",
        shear_center_ref="ex1_shear_center_closed/Shear Center Ref.xlsx"
    )
    solver = SingleClosedSectionSolver(shape)
    solver.solve()

    S_z = dimensions['S_z']
    assert expand(nsimplify(shape.ey, rational=True) - 20/7*a) == 0



def test_ex2_shear_center_closed():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': a, 'h': a})
    shape = Shape(
        nodes="ex2_shear_center_closed/Nodes.xlsx",
        elements="ex2_shear_center_closed/Elements.xlsx",
        forces="ex2_shear_center_closed/Forces.xlsx",
        sections="ex2_shear_center_closed/ClosedSections.xlsx",
        shear_center_ref="ex2_shear_center_closed/Shear Center Ref.xlsx"
    )
    solver = SingleClosedSectionSolver(shape)
    solver.solve()

    S_y = dimensions['S_y']
    assert expand(nsimplify(shape.ez, rational=True) - 22/63*a) == 0
