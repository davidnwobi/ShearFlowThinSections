from data_setup import dimensions
from open_section_solver import OpenSectionSolver
from sympy import integrate, simplify, sqrt


class SingleClosedSectionSolver(OpenSectionSolver):

    def __init__(self, shape):
        super().__init__(shape)

    def solve(self):
        self.solve_open_section()
        self.calculate_qo_moment_balance()
        self.solve_for_shear_center()
        self.update_solution()
        return self.shape



    def calculate_qo_moment_balance(self):
        moments = self.calculate_moments()
        print(f'Moments = {moments}')
        qo = - moments / (2 * self.shape.areas[0])
        print(f'qo = {qo}')
        self.shape.qo = qo

    def calculate_qo_rate_of_twist(self):
        Q = 0
        s = 0
        for element in self.shape.elements.values():
            Q += element.Q
            s += element.length
        qo = -Q / s
        print(f'qo = {qo}')
        return qo

    def update_solution(self):
        for element in self.shape.elements.values():
            if element.pos not in self.shape.sections[0]:
                continue
            element.qs += self.shape.qo
            element.Q = integrate(element.qs, (element.S, 0, element.length))
            element.Tau = element.qs / element.t
            element.qs = simplify(element.qs)
            element.Q = simplify(element.Q)
            element.Tau = simplify(element.Tau)

