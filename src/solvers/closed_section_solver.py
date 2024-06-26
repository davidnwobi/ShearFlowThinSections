from ..utils.data_setup import dimensions
from .open_section_solver import OpenSectionSolver
from sympy import integrate, simplify


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
        moments = self.calculate_moments_about_shear_force_centre()
        print(f'Moments = {moments}')
        qo = - moments / (2 * self.shape.areas[0])
        print(f'qo = {qo}')
        self.shape.qo = qo

    def calculate_moments_about_shear_force_centre(self):
        return self.calculate_moments_about_ref_point(self.shape.force_dict['S_z'][1], self.shape.force_dict['S_y'][1])

    def calculate_qo_rate_of_twist(self):
        Q = 0
        s = 0
        for element in self.shape.elements.values():
            Q += element.Q*dimensions['t']/element.t
            s += element.length*dimensions['t']/element.t
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

    def solve_for_shear_center(self):
        qo = self.calculate_qo_rate_of_twist()
        moment = self.calculate_moments_about_shear_center_ref()
        print(f'qo solve_for_shear_center= {qo}')
        print(f'Moment solve_for_shear_center= {moment}')
        print(f'Area = {self.shape.areas[0]}')
        if dimensions['S_z'] != 0:
            self.shape.ey = (moment + 2 * self.shape.areas[0] * qo) / dimensions['S_z']
        if dimensions['S_y'] != 0:
            self.shape.ez = -(moment + 2 * self.shape.areas[0] * qo) / dimensions['S_y']
        self.shape.ey = simplify(self.shape.ey)
        self.shape.ez = simplify(self.shape.ez)
