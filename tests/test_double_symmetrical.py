
from src import Shape
from sympy import symbols, expand
from src import dimensions
from src import SectionType, SectionSolver


def test_ex1_double_symmetrical():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 4 * a, 'h': a})
    shape = Shape(nodes='ex1_double_symmetrical/Nodes.xlsx',
                    elements='ex1_double_symmetrical/Elements.xlsx',
                    forces='ex1_double_symmetrical/Forces.xlsx',
                    sections='ex1_double_symmetrical/ClosedSections.xlsx',
                    shear_center_ref='ex1_double_symmetrical/Shear Center Ref.xlsx')
    SectionSolver.solve(shape, SectionType.DOUBLE_SYMMETRIC)

    S_z = dimensions['S_z']
    assert expand(shape.qo[0] - 59*S_z/(108*a)) == 0
    assert expand(shape.qo[1] - (-5*S_z/(108*a))) == 0