from data_setup import dimensions
from open_section_solver import OpenSectionSolver
from sympy import integrate, simplify, sqrt, Eq, solve, nsimplify, sign, Abs
from sympy.vector import CoordSys3D
from sympy import symbols


class DoubleSymmetricClosedSectionSolver(OpenSectionSolver):

    def __init__(self, shape):
        super().__init__(shape)

    def solve(self):
        self.solve_open_section()
        self.solve_for_qo()
        self.update_solution()
        return self.shape

    def calculate_rate_of_twist_eqns(self):
        for i, section in enumerate(self.shape.sections):
            section_rate_of_twist_eqn = 0
            self.shape.rate_of_twist_eqn.append(self.create_rate_of_twist_eqn(i, section))

    def create_rate_of_twist_eqn(self, section_no, section):
        section_rate_of_twist_eqn = 0
        qI, qII = symbols('qI qII', real=True)
        for element in self.shape.elements.values():
            if element.pos not in section:
                continue
            section_rate_of_twist_eqn += (Abs(element.Q) / element.t *
                                          sign(self.calculate_element_moment_about_ref_point(element,
                                                                                             self.shape.section_Yc[section_no],
                                                                                             self.shape.section_Zc[section_no])))
            if section_no == 0 and element.pos in self.shape.sections[0]:
                section_rate_of_twist_eqn += qI / element.t * element.length
            elif section_no == 1 and element.pos in self.shape.sections[0]:
                section_rate_of_twist_eqn -= qI / element.t * element.length

            if section_no == 1 and element.pos in self.shape.sections[1]:
                section_rate_of_twist_eqn += qII / element.t * element.length
            elif section_no == 0 and element.pos in self.shape.sections[1]:
                section_rate_of_twist_eqn -= qII / element.t * element.length
        return section_rate_of_twist_eqn

    def update_solution(self):
        for i, section in enumerate(self.shape.sections):
            for element in self.shape.elements.values():
                if element.pos not in section:
                    continue
                moment_sign = sign(self.calculate_element_moment_about_ref_point(element,
                                                                                 self.shape.section_Yc[
                                                                                     i],
                                                                                 self.shape.section_Zc[

                                                                                    i], Qb=True))
                element.qs += self.shape.qo[i]*moment_sign*sign(element.Qb)
                element.Q = integrate(element.qs, (element.S, 0, element.length))
                element.Tau = element.qs / element.t
                element.qs = simplify(element.qs)
                element.Q = simplify(element.Q)
                element.Tau = simplify(element.Tau)

    def calculate_moments_about_ref_point(self, y_1, z_1):
        # Calculate the moments
        M = 0
        for i, section in enumerate(self.shape.sections):
            for index, element in self.shape.elements.items():
                if element.pos not in section:
                    continue
                M += self.calculate_element_moment_about_ref_point(element, y_1, z_1)
        M = simplify(M)
        return M

    @staticmethod
    def calculate_element_moment_about_ref_point(element, y_1, z_1, Qb=False):
        M = 0
        N = CoordSys3D('N')
        y_2 = element.Node1.y + (element.Node2.y - element.Node1.y) / 2
        z_2 = element.Node1.z + (element.Node2.z - element.Node1.z) / 2
        F = 0
        if Qb:
            F = element.Qb * (element.cos() * N.j + element.sin() * N.k)
        else:
            F = element.Q * (element.cos() * N.j + element.sin() * N.k)
        d = (y_2 - y_1) * N.j + (z_2 - z_1) * N.k
        M -= (F.cross(d)).dot(N.i)
        return M

    def solve_for_qo(self):
        self.calculate_rate_of_twist_eqns()
        eq1 = Eq(self.shape.rate_of_twist_eqn[0] - self.shape.rate_of_twist_eqn[1], 0)
        eq1 = nsimplify(eq1)
        eq1 = simplify(eq1)

        qI, qII = symbols('qI qII', real=True)
        moment = 0
        result = None

        if dimensions['S_z'] != 0:
            moment = self.calculate_moments_about_ref_point(0, 0)
            eq2 = Eq(dimensions['S_z'] * self.shape.force_dict['S_z'][1] - (
                        moment + 2 * (self.shape.areas[0] * qI + self.shape.areas[1] * qII)), 0)
            eq2 = simplify(eq2)
            result = solve((eq1, eq2), (qI, qII))

        if dimensions['S_y'] != 0:
            moment = self.calculate_moments_about_ref_point(0, 0)
            eq2 = Eq(-dimensions['S_y'] * self.shape.force_dict['S_y'][1],
                     moment + 2 * (self.shape.areas[0] * qI + self.shape.areas[1] * qII))
            eq2 = simplify(eq2)

            result = solve([eq1, eq2], (qI, qII))

        qI = result[qI]
        qII = result[qII]
        qI = nsimplify(qI)
        qII = nsimplify(qII)
        qI = simplify(qI)
        qII = simplify(qII)

        print('qI =', qI)
        print('qII =', qII)
        self.shape.qo = [qI, qII]
