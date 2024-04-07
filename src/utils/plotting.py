from .data_setup import dimensions
from sympy import symbols


def draw_shape(shape):
    import matplotlib.pyplot as plt
    import numpy as np
    b = symbols('b', real=True, positive=True)
    fig, ax = plt.subplots()
    for element in shape.elements.values():
        print(element.Node1.y.as_coefficients_dict())
        x1 = sum(list(element.Node1.y.as_coefficients_dict().values()))
        y1 = sum(list(element.Node1.z.as_coefficients_dict().values()))
        x2 = sum(list(element.Node2.y.as_coefficients_dict().values()))
        y2 = sum(list(element.Node2.z.as_coefficients_dict().values()))
        ax.plot([x1, x2], [y1, y2], 'k-', lw=element.t.subs(dimensions['t'], 1))

    ax.axis('equal')
    ax.invert_yaxis()
    ax.invert_xaxis()
    plt.show()
