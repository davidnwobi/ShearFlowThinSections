from shape import Shape
from sympy import symbols, expand
from data_setup import dimensions
from open_section_solver import OpenSectionSolver


def test_ex1_shear_center_open():
    b, h = symbols('b h', real=True, positive=True)
    dimensions.update({'b': b, 'h': h})
    shape = Shape()
    solver = OpenSectionSolver(shape)
    solver.solve()

    S_z = dimensions['S_z']
    assert expand(shape.ey - 3 * b ** 2 / (6 * b + h)) == 0
