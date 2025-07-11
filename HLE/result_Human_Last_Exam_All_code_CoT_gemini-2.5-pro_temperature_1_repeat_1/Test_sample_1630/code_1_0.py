import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

def solve_task():
    """
    This script defines two cubic polynomials f(x) and g(x) with positive derivatives
    and finds the fixed points of their composition h(x) = f(g(x)).
    """

    # Define the polynomials f(y) and g(x)
    # g(x) = (0.01/3)x^3 - 0.2x^2 + 4.1x
    # f(y) = (0.01/3)y^3 + 0.5y
    def g(x):
        return (0.01/3.0)*x**3 - 0.2*x**2 + 4.1*x

    def f(y):
        return (0.01/3.0)*y**3 + 0.5*y

    # Their derivatives must be positive for all x
    # g'(x) = 0.01x^2 - 0.4x + 4.1. The discriminant is (-0.4)^2 - 4*0.01*4.1 = -0.004 < 0,
    # so g'(x) is always positive.
    # f'(y) = 0.01y^2 + 0.5, which is always positive.
    
    # The function for which we want to find the roots (h(x) - x = 0)
    def h_minus_x(x):
        return f(g(x)) - x

    # We will find the roots by searching for sign changes in a given range.
    # Plotting the function helps visualize the roots and choose a search interval.
    x_vals = np.linspace(-15, 45, 4000)
    y_vals = h_minus_x(x_vals)

    # Find brackets for the roots (where the function crosses zero)
    roots = []
    for i in range(len(x_vals) - 1):
        if np.sign(y_vals[i]) != np.sign(y_vals[i+1]):
            # Use a numerical solver to find the precise root in the bracket
            try:
                root = brentq(h_minus_x, x_vals[i], x_vals[i+1])
                # Check if root is not too close to the last one found
                if not any(np.isclose(root, r) for r in roots):
                    roots.append(root)
            except ValueError:
                # brentq fails if y_vals[i] or y_vals[i+1] is exactly zero
                pass

    print(f"Found {len(roots)} fixed points for f(g(x)) = x.")
    print("The fixed points are:")
    for root in sorted(roots):
        # h_val should be very close to root_val
        h_val = f(g(root))
        print(f"f(g({root:.9f})) = {h_val:.9f}")

solve_task()
