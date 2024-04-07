from dataclasses import dataclass
from sympy import Mul


@dataclass
class Node:
    pos: int
    y: Mul
    z: Mul

