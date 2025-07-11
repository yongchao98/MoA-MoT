import numpy as np
from scipy.integrate import quad

def solve_and_print():
    """
    Solves the problem by finding the circle's equation, calculating alpha,
    and then determining the probability P(X < 3).
    """
    # 1. Define the parameters of the circle found from the points A, B, and C.
    # Center (h, k), and radius squared r_sq.
    h = 5.5
    k = -1.5
    r_sq = 22.5
    
    # Print the equation of the circle, as requested.
    print(f"The equation of the circle is (x - {h})^2 + (y - ({k}))^2 = {r_sq}")

    # 2. Define the function f(x) representing the upper arc of the circle.
    def f(x):
        # f(x) = k + sqrt(r^2 - (x-h)^2) - but our k is negative.
        # (y-k)^2 -> y- (-1.5) = y+1.5
        # y = -1.5 + sqrt(...)
        return -1.5 + np.sqrt(r_sq - (x - h)**2)

    # 3. Calculate alpha.
    # The total area under f(x) from x=1 to x=10.
    # The result of quad is a tuple (value, error). We only need the value.
    total_integral, _ = quad(f, 1, 10)
    
    # alpha is the normalization constant.
    alpha = 1.0 / total_integral

    # 4. Calculate the probability P(X < 3).
    # This is the area under the PDF from x=1 to x=3.
    # The random variable X is defined on [1, 10], so P(X < 3) is the integral from 1 to 3.
    partial_integral, _ = quad(f, 1, 3)
    
    # P(X < 3) = alpha * partial_integral
    probability = partial_integral / total_integral

    print(f"The value of alpha is: {alpha:.6f}")
    print(f"The value of P(X < 3) is: {probability:.6f}")

solve_and_print()
<<<0.028692, 0.147122>>>