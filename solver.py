from abc import ABC, abstractmethod


class Solver(ABC):
    def __init__(self, shape):
        self.shape = shape

    @abstractmethod
    def solve(self):
        pass
