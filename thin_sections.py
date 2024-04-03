import pandas as pd

from shape import Shape
from sympy import *
from data_setup import dimensions
from open_section_solver import OpenSectionSolver
from closed_section_solver import SingleClosedSectionSolver

a = symbols('a', real=True, positive=True)
dimensions.update({'b': 2 * a, 'h': 2 * a})
shape = Shape()
solver = SingleClosedSectionSolver(shape)
solver.solve()

print('Shape Properties:')
print(f'dy = {shape.dy}')
print(f'dz = {shape.dz}')
print(f'Iy = {shape.Iy}')
print(f'Iz = {shape.Iz}')
print(f'Iyz = {shape.Iyz}')
print(f'Ey = {shape.ey}')
print(f'Ez = {shape.ez}\n')
print('Element Properties:')
solution = pd.DataFrame(columns=['Element', 'l', 't', 'ty', 'tz', 'Int_ty', 'Int_tz', 'qs', 'Q', 'Tau'])
for element in shape.elements.values():
    mini_df = pd.DataFrame(
        {'Element': [element.pos], 'l': [element.length], 't': [element.t], 'ty': [element.ty], 'tz': [element.tz],
         'Int_ty': [element.int_ty], 'Int_tz': [element.int_tz], 'qs': [element.qs], 'Q': [element.Q],
         'Tau': [element.Tau]})
    solution = pd.concat([solution, mini_df])

    print(f'Element {element.pos}:')
    print(f'l = {element.length}')
    print(f't = {element.t}')
    print(f'ty = {element.ty}')
    print(f'tz = {element.tz}')
    print(f'Int_ty = {element.int_ty}')
    print(f'Int_tz = {element.int_tz}')
    print(f'Int_ty at Node1 = {element.int_ty.subs(element.S, 0)}')
    print(f'Int_ty at Node2 = {element.int_ty.subs(element.S, element.length)}')
    print(f'Int_tz at Node1 = {element.int_tz.subs(element.S, 0)}')
    print(f'Int_tz at Node2 = {element.int_tz.subs(element.S, element.length)}')
    print(f'qb = {element.qb}')
    print(f'qb at Node1 = {element.qb.subs(element.S, 0)}')
    print(f'qb at Node2 = {element.qb.subs(element.S, element.length)}')
    print(f'Qb = {element.Qb}')
    print(f'qs = {element.qs}')
    print(f'qs at Node1 = {element.qs.subs(element.S, 0)}')
    print(f'qs at Node2 = {element.qs.subs(element.S, element.length)}')
    print(f'Q = {element.Q}')
    print(f'Tau = {element.Tau}\n')


def obtain_max_value_of_equation(equation, variable):
    derivative = diff(equation, variable)
    if derivative == 0:
        return (equation, None)
    critical_points = solve(derivative, variable)
    max_value = equation.subs(variable, 0)
    max_value_pos = 0
    if len(critical_points) > 0:
        for point in critical_points:
            value = equation.subs(variable, point)
            if abs(value) > abs(max_value):
                max_value = value
                max_value_pos = point
    value = equation.subs(variable, element.length)
    if abs(value) > abs(max_value):
        max_value = value
        max_value_pos = element.length
    return nsimplify(max_value, rational=True, tolerance=1e-14), nsimplify(max_value_pos, rational=True, tolerance=1e-14)


for element in shape.elements.values():
    print(f'Maximum Value of qs for Element {element.pos}:')
    max_qs = obtain_max_value_of_equation(element.qs, element.S)
    print(f'Max_qs = {max_qs[0]} at S = {max_qs[1]} \n')