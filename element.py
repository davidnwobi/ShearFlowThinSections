from sympy import Mul, sqrt, symbols, Symbol
from node import Node

class Element:
    def __init__(self, pos: int,  Node1: 'Node', Node2: 'Node', t: Symbol):
        self.pos = pos
        self.Node1 = Node1
        self.Node2 = Node2
        self.t = t
        self.S = symbols(f'S_{Node1.pos}{Node2.pos}')
        self.y = None
        self.z = None
        self.ty = None
        self.tz = None
        self.int_ty = None
        self.int_tz = None
        self.qb = None
        self.qs = None
        self.Qb = None
        self.Q = None
        self.Taub = None
        self.Tau = None

    @property
    def length(self) -> Mul:
        return sqrt((self.Node1.y - self.Node2.y) ** 2 + (self.Node1.z - self.Node2.z) ** 2)

    def __repr__(self):
        return f'Element {self.Node1.pos} -> {self.Node2.pos}, t={self.t}'

    def sin(self):
        return (self.Node2.z - self.Node1.z) / self.length

    def cos(self):
        return (self.Node2.y - self.Node1.y) / self.length

    def tan(self):
        return (self.Node2.z - self.Node1.z) / (self.Node2.y - self.Node1.y)


