from shape import Shape
from sympy import symbols, expand
from data_setup import dimensions
from open_section_solver import OpenSectionSolver
from plotting import draw_shape


def test_ex1_shear_center_open():
    b, h = symbols('b h', real=True, positive=True)
    dimensions.update({'b': b, 'h': h})
    shape = Shape(
        nodes='ex1_shear_center_open/Nodes.xlsx',
        elements='ex1_shear_center_open/Elements.xlsx',
        forces='ex1_shear_center_open/Forces.xlsx',
        sections='ex1_shear_center_open/ClosedSections.xlsx',
        shear_center_ref='ex1_shear_center_open/Shear Center Ref.xlsx')
    solver = OpenSectionSolver(shape)
    solver.solve()

    S_z = dimensions['S_z']
    assert expand(shape.ey - 3 * b ** 2 / (6 * b + h)) == 0
