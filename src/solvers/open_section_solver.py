from sympy import integrate, simplify, nsimplify
from src.utils.data_setup import dimensions
from .solver import Solver
from sympy.vector import CoordSys3D


class OpenSectionSolver(Solver):
    def __init__(self, shape):
        super().__init__(shape)

    def solve(self):
        self.solve_open_section()
        self.solve_for_shear_center()
        return self.shape

    # TODO: Properly evaluate the moments

    def solve_open_section(self):
        can_continue = True
        while can_continue:
            solved_on_this_iteration = 0
            for index, element in self.shape.elements.items():
                if not self.element_is_solvable(element):
                    continue
                # Integrate the element's ty and tz
                self.solve_for_integrals(element)
                self.solve_for_shear_flow(element)
                solved_on_this_iteration += 1
            if solved_on_this_iteration == 0:
                can_continue = False

    def calculate_moments_about_shear_center_ref(self):
        return self.calculate_moments_about_ref_point(self.shape.shear_center_y, self.shape.shear_center_z)

    def calculate_moments_about_ref_point(self, y_1, z_1):
        # Calculate the moments
        M = 0
        for i, section in enumerate(self.shape.sections):
            for index, element in self.shape.elements.items():
                if element.pos not in section:
                    continue
                M += self.calculate_element_moment_about_ref_point(element, y_1, z_1)
        M = nsimplify(M)
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
        M -= nsimplify((F.cross(d)).dot(N.i))
        return M

    def solve_for_shear_center(self):
        moment = self.calculate_moments_about_shear_center_ref()

        if dimensions['S_z'] != 0:
            self.shape.ey = moment / dimensions['S_z']
        if dimensions['S_y'] != 0:
            self.shape.ez = - moment / dimensions['S_y']
        self.shape.ey = simplify(self.shape.ey)
        self.shape.ez = simplify(self.shape.ez)

    def solve_for_integrals(self, element):
        element.y = element.Node1.y + (element.Node2.y - element.Node1.y) / element.length * element.S
        element.z = element.Node1.z + (element.Node2.z - element.Node1.z) / element.length * element.S
        element.y = simplify(element.y)
        element.z = simplify(element.z)

        element.ty = element.t * element.y
        element.tz = element.t * element.z
        element.ty = simplify(element.ty)
        element.tz = simplify(element.tz)

        element.int_ty = integrate(element.ty, element.S)
        element.int_tz = integrate(element.tz, element.S)
        element.int_ty = simplify(element.int_ty)
        element.int_tz = simplify(element.int_tz)

        self.solve_for_constants(element)

        element.int_ty = nsimplify(element.int_ty, tolerance=1e-10, rational=True)
        element.int_tz = nsimplify(element.int_tz, tolerance=1e-10, rational=True)
        element.int_ty = simplify(element.int_ty)
        element.int_tz = simplify(element.int_tz)

    def solve_for_shear_flow(self, element):
        pt1_num = self.shape.Iz * dimensions['S_z'] - self.shape.Iyz * dimensions['S_y']
        pt1_denom = self.shape.Iy * self.shape.Iz - self.shape.Iyz ** 2
        pt2_num = self.shape.Iy * dimensions['S_y'] - self.shape.Iyz * dimensions['S_z']
        pt2_denom = self.shape.Iy * self.shape.Iz - self.shape.Iyz ** 2

        if pt1_denom == 0 or pt2_denom == 0:
            ValueError('Something went wrong. A denominator are zero.')

        if pt1_num == 0:
            element.qs = -(pt2_num / pt2_denom) * element.int_ty
        elif pt2_num == 0:
            element.qs = -(pt1_num / pt1_denom) * element.int_tz
        else:
            element.qs = -(pt1_num / pt1_denom * element.int_tz) - (pt2_num / pt2_denom * element.int_ty)

        element.Q = integrate(element.qs, (element.S, 0, element.length))
        element.Tau = element.qs / element.t

        element.qs = nsimplify(element.qs, tolerance=1e-10, rational=True)
        element.qb = element.qs
        element.Q = nsimplify(element.Q, tolerance=1e-10, rational=True)
        element.Qb = element.Q
        element.Tau = nsimplify(element.Tau, tolerance=1e-10, rational=True)
        element.Taub = element.Tau

        element.qs = simplify(element.qs)
        element.qb = simplify(element.qb)
        element.Q = simplify(element.Q)
        element.Qb = simplify(element.Qb)
        element.Tau = simplify(element.Tau)
        element.Taub = simplify(element.Taub)

    def solve_for_constants(self, element):
        # Locate which of the two nodes has a boundary condition
        solvable_node = self.locate_solvable_node(element)

        int_ty_at_node = 0
        int_tz_at_node = 0

        # if the element is an end element and the solvable node is not a terminal node

        # TODO: Work on this condition.It seems redundant. It assumes that shear flow always flows in from the free end
        if len(self.shape.node_graph[solvable_node]) == 1:  # Check if the node is a terminal node
            if solvable_node != min(element.Node1.pos, element.Node2.pos):  # Check flow direction
                element.int_ty += 0 - element.int_ty.subs(element.S, element.length)
                element.int_tz += 0 - element.int_tz.subs(element.S, element.length)
            return True
        else:
            lhs_ty = 0
            rhs_ty = 0
            lhs_tz = 0
            rhs_tz = 0

            # Obtain the sum of the integrals of ty and tz going into and out of the solvable node
            for node in self.shape.node_graph[solvable_node]:
                if node < solvable_node:
                    element = self.shape.elements[(node, solvable_node)]
                    lhs_ty += element.int_ty.subs(element.S, element.length)
                    lhs_tz += element.int_tz.subs(element.S, element.length)
                else:
                    element = self.shape.elements[(solvable_node, node)]
                    rhs_ty += element.int_ty.subs(element.S, 0)
                    rhs_tz += element.int_tz.subs(element.S, 0)
            # Solve the eqn X_before - X_after = int_X
            if solvable_node == min(element.Node1.pos, element.Node2.pos):
                int_ty_at_node = lhs_ty - rhs_ty
                int_tz_at_node = lhs_tz - rhs_tz
            else:
                # Solve the eqn X_after - X_before = int_X
                int_ty_at_node = rhs_ty - lhs_ty
                int_tz_at_node = rhs_tz - lhs_tz

            # Apply the boundary condition to the solvable node
            if solvable_node == min(element.Node1.pos, element.Node2.pos):
                element.int_ty += int_ty_at_node
                element.int_tz += int_tz_at_node
                return True
            else:
                int_ty_at_node += int_ty_at_node - element.int_ty.subs(element.S, element.length)
                int_tz_at_node += int_tz_at_node - element.int_tz.subs(element.S, element.length)
                return False

    def element_is_solvable(self, element):
        if element.int_ty is not None or element.int_tz is not None:
            return False
        if len(self.shape.node_graph[element.Node1.pos]) == 1 or len(self.shape.node_graph[element.Node2.pos]) == 1:
            return True

        solvable = []
        for element_node in [element.Node1.pos, element.Node2.pos]:
            for other_node in self.shape.node_graph[element_node]:
                i = min(element_node, other_node)
                j = max(element_node, other_node)
                if (i, j) != (element.Node1.pos, element.Node2.pos) and (
                        self.shape.elements[(i, j)].int_ty is None or self.shape.elements[(i, j)].int_tz is None):
                    solvable.append(False)
                    break
            else:
                solvable.append(True)
        return any(solvable)

    # TODO: Keep track of solvable nodes rather than checking all nodes
    def locate_solvable_node(self, element):
        if len(self.shape.node_graph[element.Node1.pos]) == 1:
            return element.Node1.pos
        if len(self.shape.node_graph[element.Node2.pos]) == 1:
            return element.Node2.pos
        for element_node in [element.Node1.pos, element.Node2.pos]:
            for other_node in self.shape.node_graph[element_node]:
                i = min(element_node, other_node)
                j = max(element_node, other_node)
                if self.shape.elements[(i, j)].int_ty is None or self.shape.elements[(i, j)].int_tz is None:
                    break
            else:
                return element_node
