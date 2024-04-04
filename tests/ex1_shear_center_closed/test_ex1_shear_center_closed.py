from shape import Shape
from sympy import symbols, expand, nsimplify
from data_setup import dimensions
from closed_section_solver import SingleClosedSectionSolver


def test_ex1_shear_center_closed():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 5 * a, 'h': 4.8 * a})
    shape = Shape()
    solver = SingleClosedSectionSolver(shape)
    solver.solve()

    S_z = dimensions['S_z']
    assert expand(nsimplify(shape.ey, rational=True) - 20/7*a) == 0



