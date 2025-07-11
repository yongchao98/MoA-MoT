import math
from scipy.integrate import quad

def solve_and_print():
    """
    Calculates the values of alpha and P(X<3) based on the problem description.
    """
    # Parameters of the circle equation derived from points A(1,0), B(10,0), and C(4,3)
    # Equation: (x-h)^2 + (y-k)^2 = r2
    h = 5.5
    k = -1.5
    r2 = 22.5
    
    print("Step 1: Determine the equation of the circle and the function f(x).")
    print(f"The equation of the circle is (x - {h})^2 + (y - ({k}))^2 = {r2}")
    print(f"The function y = f(x) for the upper arc is f(x) = {k} + sqrt({r2} - (x - {h})^2)")
    print("-" * 20)

    # Define the function f(x) which describes the upper arc of the circle.
    def f(x):
        # radicand must be non-negative.
        radicand = r2 - (x - h)**2
        # The integration limits are within the domain of the function, so this check is a safeguard.
        if radicand < 0:
            return 0.0
        return k + math.sqrt(radicand)

    # Calculate the total integral of f(x) over the domain [1, 10]
    # This integral corresponds to the area under the curve f(x)
    total_integral, total_err = quad(f, 1, 10)

    # Calculate alpha such that the total probability is 1
    # alpha = 1 / total_integral
    alpha = 1.0 / total_integral

    # Calculate the partial integral for P(X<3), which means integrating from 1 to 3
    partial_integral, partial_err = quad(f, 1, 3)

    # Calculate the probability P(X<3)
    # P(X<3) = alpha * partial_integral
    p_x_lt_3 = partial_integral / total_integral
    
    print("Step 2: Calculate alpha and P(X < 3).")
    print(f"The value of alpha is: {alpha}")
    print(f"The value of P(X < 3) is: {p_x_lt_3}")

solve_and_print()