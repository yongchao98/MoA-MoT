import numpy as np
from scipy.integrate import quad

def solve_circle_pdf():
    """
    Solves for the parameters of the PDF and calculates the required probability.
    """
    # Step 1 & 2: Determine the circle's equation and the function f(x)
    # The center of the circle (h,k) and the square of its radius (r_squared)
    # have been derived from the points A(1,0), B(10,0), and C(4,3).
    h = 5.5
    k = -1.5
    r_squared = 22.5

    # The function f(x) corresponds to the upper arc of the circle.
    # y = f(x) = k + sqrt(r^2 - (x-h)^2)
    def f(x):
        radicand = r_squared - (x - h)**2
        # Handle potential floating-point inaccuracies at the domain boundaries
        if isinstance(x, (np.ndarray, list, tuple)):
            radicand[radicand < 0] = 0
        elif radicand < 0:
            radicand = 0
        return k + np.sqrt(radicand)

    # Step 3: Calculate the normalization constant alpha
    # The integral of the PDF d_X(x) = alpha * f(x) over [1, 10] must be 1.
    # So, alpha = 1 / integral(f(x) from 1 to 10)
    total_area, _ = quad(f, 1, 10)
    alpha = 1.0 / total_area

    # Step 4: Calculate the probability P(X < 3)
    # P(X < 3) = integral(alpha * f(x) from 1 to 3)
    area_less_than_3, _ = quad(f, 1, 3)
    probability = alpha * area_less_than_3

    # Output the results
    print("--- Circle Equation Parameters ---")
    print(f"The equation of the circle is (x - h)^2 + (y - k)^2 = r^2")
    print(f"h = {h}")
    print(f"k = {k}")
    print(f"r^2 = {r_squared}")
    print("\n--- Final Results ---")
    print(f"The value of alpha is: {alpha}")
    print(f"The value of P(X < 3) is: {probability}")
    
    # Store results for final answer block
    global alpha_val, prob_val
    alpha_val = alpha
    prob_val = probability

solve_circle_pdf()

<<<0.04683108615481745, 0.07658937013149557>>>