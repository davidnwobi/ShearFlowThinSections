from shape import Shape
from sympy import symbols, expand, nsimplify
from data_setup import dimensions
from closed_section_solver import SingleClosedSectionSolver


def test_ex1_shear_center_closed():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': a, 'h': a})
    shape = Shape()
    solver = SingleClosedSectionSolver(shape)
    solver.solve()

    S_y = dimensions['S_y']
    assert expand(nsimplify(shape.ez, rational=True) - 22/63*a) == 0



