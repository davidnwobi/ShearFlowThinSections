from ..data.shape import Shape
from ..solvers.open_section_solver import OpenSectionSolver
from ..solvers.closed_section_solver import SingleClosedSectionSolver
from ..solvers.double_symmetric_closed_section_solver import DoubleSymmetricClosedSectionSolver
from enum import Enum


class SectionType(str, Enum):
    OPEN = 'Open'
    SINGLE_CLOSED = 'Single Closed'
    DOUBLE_SYMMETRIC = 'Double Symmetric'


class SectionSolver:
    @staticmethod
    def solve(shape: Shape, section_type: SectionType) -> None:
        if section_type == SectionType.OPEN:
            solver = OpenSectionSolver(shape)
        elif section_type == SectionType.SINGLE_CLOSED:
            solver = SingleClosedSectionSolver(shape)
        elif section_type == SectionType.DOUBLE_SYMMETRIC:
            solver = DoubleSymmetricClosedSectionSolver(shape)
        else:
            raise ValueError('Invalid Section Type')
        solver.solve()
