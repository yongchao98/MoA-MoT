import numpy as np
from scipy.integrate import quad

def solve_problem():
    """
    This function solves the problem by finding the values of alpha and P(X < 3).
    """

    # The equation of the circle passing through A(1,0), B(10,0), and C(4,3) is
    # (x - h)^2 + (y - k)^2 = r_squared.
    # We determined the center (h, k) and the radius squared r_squared.
    h = 5.5
    k = -1.5
    r_squared = 22.5

    # The function y = f(x) represents the upper arc of the circle.
    # The equation for y is derived as y = k + sqrt(r_squared - (x - h)^2).
    # Since we have (y-k)^2, and y+1.5=sqrt(), it's y=-1.5+sqrt().
    # So our formula should be y = -1.5 + sqrt(...) which corresponds to the value of k.
    print(f"The equation for the upper arc of the circle is y = f(x) = {k} + sqrt({r_squared} - (x - {h})^2)")

    def f(x):
        """The function defining the upper arc of the circle."""
        radicand = r_squared - (x - h)**2
        # Ensure the value inside the square root is non-negative to avoid domain errors
        # with floating point inaccuracies near the boundaries.
        if isinstance(radicand, (np.ndarray, list, tuple)):
            radicand = np.maximum(radicand, 0)
        elif radicand < 0:
            radicand = 0
        return k + np.sqrt(radicand)

    # To find alpha, we need to ensure that the integral of the PDF over [1, 10] is 1.
    # alpha * integral(f(x) dx) from 1 to 10 = 1
    # So, alpha = 1 / integral(f(x) dx).
    total_integral, _ = quad(f, 1, 10)

    if total_integral <= 0:
        print("Error: The integral of f(x) is not positive. The PDF is not valid.")
        return

    alpha = 1 / total_integral

    # To find P(X < 3), we integrate the PDF from 1 to 3.
    # P(X < 3) = integral(alpha * f(x) dx) from 1 to 3
    # P(X < 3) = alpha * integral(f(x) dx) from 1 to 3
    partial_integral, _ = quad(f, 1, 3)

    probability = alpha * partial_integral

    print(f"\nThe calculated value of alpha is: {alpha}")
    print(f"The calculated value of P(X < 3) is: {probability}")

# Run the solver
solve_problem()