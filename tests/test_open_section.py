from src import Shape
from sympy import symbols, expand
from src import dimensions
from src import SectionType, SectionSolver


def test_ex1_open_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': a, 'h': 2 * a})
    shape = Shape(nodes='ex1_open_section/Nodes.xlsx',
                  elements='ex1_open_section/Elements.xlsx',
                  forces="ex1_open_section/Forces.xlsx",
                  sections="ex1_open_section/ClosedSections.xlsx",
                  shear_center_ref="ex1_open_section/Shear Center Ref.xlsx")
    SectionSolver.solve(shape, SectionType.OPEN)

    t = dimensions['t']
    Sz = dimensions['S_z']
    S01 = shape.elements[(0, 1)].S
    l01 = shape.elements[(0, 1)].length
    assert shape.Iy == 14 * a ** 3 * t / 3
    assert shape.elements[(0,1)].tz == -2*a*t
    assert shape.elements[(0,1)].int_tz == -2*a*t*S01
    assert shape.elements[(0,1)].qs.subs(S01, a) == 3*Sz/(7*a)
    assert shape.elements[(0,1)].Q == 3*Sz/14

    l12 = shape.elements[(1,2)].length
    S12 = shape.elements[(1, 2)].S
    assert expand(shape.elements[(1, 2)].tz - t*(S12-a)) == 0
    assert expand(shape.elements[(1, 2)].int_tz - (t*(S12**2/2-a*S12)-2*a**2*t)) == 0
    assert expand(shape.elements[(1, 2)].qs.subs(S12, a)-15*Sz/(28*a)) == 0
    assert expand(shape.elements[(1, 2)].Q - Sz) == 0

    l23 = shape.elements[(2, 3)].length
    S23 = shape.elements[(2, 3)].S
    assert shape.elements[(2, 3)].tz == 2*a*t
    assert expand(shape.elements[(2, 3)].int_tz - 2*a*t*(S23-a)) == 0
    assert expand(shape.elements[(2, 3)].qs.subs(S23, 0) - 3*Sz/(7*a)) == 0
    assert expand(shape.elements[(2, 3)].Q - 3*Sz/14) == 0


def test_ex2_open_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': a, 'h': 2 * a})
    shape = Shape(nodes='ex2_open_section/Nodes.xlsx',
                    elements='ex2_open_section/Elements.xlsx',
                    forces="ex2_open_section/Forces.xlsx",
                    sections="ex2_open_section/ClosedSections.xlsx",
                    shear_center_ref="ex2_open_section/Shear Center Ref.xlsx")
    SectionSolver.solve(shape, SectionType.OPEN)
    t = dimensions['t']
    Sz = dimensions['S_z']
    S01 = shape.elements[(0, 1)].S
    l01 = shape.elements[(0, 1)].length
    assert shape.Iy == 14 * a ** 3 * t / 3

    assert expand(shape.elements[(0, 1)].tz - (-2 * a * t)) == 0
    assert expand(shape.elements[(0, 1)].int_tz - (-2 * a * t * S01)) == 0
    assert shape.elements[(0,1)].qs.subs(S01, a) == 3*Sz/(7*a)
    assert shape.elements[(0,1)].Q == 3*Sz/14

    S13 = shape.elements[(1, 3)].S
    l13 = shape.elements[(1, 3)].length

    assert expand(shape.elements[(1, 3)].tz - t * (S13 - a)) == 0
    assert expand(shape.elements[(1, 3)].int_tz - (t * (S13 ** 2 / 2 - a * S13) - 2 * a ** 2 * t)) == 0
    assert expand(shape.elements[(1, 3)].qs.subs(S13, a) - 15 * Sz / (28 * a)) == 0
    assert expand(shape.elements[(1, 3)].Q - Sz/2) == 0

    S23 = shape.elements[(2, 3)].S
    l23 = shape.elements[(2, 3)].length

    assert expand(shape.elements[(2, 3)].tz - 0) == 0
    assert expand(shape.elements[(2, 3)].int_tz - 0) == 0

    S34 = shape.elements[(3, 4)].S
    l34 = shape.elements[(3, 4)].length

    assert expand(shape.elements[(3, 4)].tz - S34 * t) == 0
    assert expand(shape.elements[(3, 4)].int_tz - (S34 ** 2 / 2 * t - 5 * a ** 2 * t / 2)) == 0
    assert expand(shape.elements[(3, 4)].qs.subs(S34, 0) - 15 * Sz / (28 * a)) == 0
    assert expand(shape.elements[(3, 4)].Q - Sz/2) == 0

    S45 = shape.elements[(4, 5)].S
    l45 = shape.elements[(4, 5)].length

    assert expand(shape.elements[(4, 5)].tz - 2 * a * t) == 0
    assert expand(shape.elements[(4, 5)].int_tz - 2 * a * t * (S45 - a)) == 0
    assert expand(shape.elements[(4, 5)].qs.subs(S45, 0) - 3*Sz/(7*a)) == 0
    assert expand(shape.elements[(4, 5)].Q - 3*Sz/14) == 0

def test_ex3_open_section():
    a = symbols('a', real=True, positive=True)
    dimensions.update({'b': 2 * a, 'h': 2 * a})
    shape = Shape(nodes='ex3_open_section/Nodes.xlsx',
                    elements='ex3_open_section/Elements.xlsx',
                    forces="ex3_open_section/Forces.xlsx",
                    sections="ex3_open_section/ClosedSections.xlsx",
                    shear_center_ref="ex3_open_section/Shear Center Ref.xlsx")
    SectionSolver.solve(shape, SectionType.OPEN)

    t = dimensions['t']
    Sz = dimensions['S_z']
    S01 = shape.elements[(0, 1)].S
    l01 = shape.elements[(0, 1)].length
    assert shape.Iy == 37 * a ** 3 * t / 6
    assert shape.Iz == 15 * a ** 3 * t / 8
    assert shape.Iyz == -11 * a ** 3 * t / 4

    assert expand(shape.elements[(0, 1)].int_ty - 4 * t * (7 / 8 * a * S01 - S01 ** 2 / 2)) == 0
    assert expand(shape.elements[(0, 1)].int_tz - (- 3 * a * t * S01)) == 0
    assert expand(shape.elements[(0, 1)].qs.subs(S01,l01) - 3 * Sz / (8 * a)) == 0

    l12 = shape.elements[(1, 2)].length
    S12 = shape.elements[(1, 2)].S

    assert expand(shape.elements[(1, 2)].int_ty - (-a / 8 * t * S12 + 3 / 2 * a ** 2 * t)) == 0
    assert expand(shape.elements[(1, 2)].int_tz - (t * (-3 * a * S12 / 4 + S12 ** 2 / 2) - 3 * a ** 2 * t)) == 0
    assert expand(shape.elements[(1, 2)].qs.subs(S12, 0) - 3 * Sz / (8 * a)) == 0
    assert expand(shape.elements[(1, 2)].qs.subs(S12, l12) - 5 * Sz / (16 * a)) == 0

    l23 = shape.elements[(2, 3)].length
    S23 = shape.elements[(2, 3)].S

    assert expand(shape.elements[(2, 3)].int_ty - (2 * t * (-a / 8 * S23 - S23 ** 2 / 2) + 5 / 4 * a ** 2 * t)) == 0
    assert expand(shape.elements[(2, 3)].int_tz - (5 / 2 * a * t * S23 - 5 / 2 * a ** 2 * t)) == 0
    assert expand(shape.elements[(2, 3)].qs.subs(S23, 0) - 5 * Sz / (16 * a)) == 0
