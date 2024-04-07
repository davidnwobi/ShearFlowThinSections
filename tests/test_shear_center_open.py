from src import Shape
from sympy import symbols, expand
from src import dimensions
from src import SectionType, SectionSolver

def test_ex1_shear_center_open():
    b, h = symbols('b h', real=True, positive=True)
    dimensions.update({'b': b, 'h': h})
    shape = Shape(
        nodes='ex1_shear_center_open/Nodes.xlsx',
        elements='ex1_shear_center_open/Elements.xlsx',
        forces='ex1_shear_center_open/Forces.xlsx',
        sections='ex1_shear_center_open/ClosedSections.xlsx',
        shear_center_ref='ex1_shear_center_open/Shear Center Ref.xlsx')
    SectionSolver.solve(shape, SectionType.OPEN)
    assert expand(shape.ey - 3 * b ** 2 / (6 * b + h)) == 0


def test_ex2_shear_center_open():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': a, 'h': 2 * a})
    shape = Shape(
        nodes='ex2_shear_center_open/Nodes.xlsx',
        elements='ex2_shear_center_open/Elements.xlsx',
        forces='ex2_shear_center_open/Forces.xlsx',
        sections='ex2_shear_center_open/ClosedSections.xlsx',
        shear_center_ref='ex2_shear_center_open/Shear Center Ref.xlsx')

    SectionSolver.solve(shape, SectionType.OPEN)
    assert expand(shape.ey.subs(a, 14) - 6) == 0
