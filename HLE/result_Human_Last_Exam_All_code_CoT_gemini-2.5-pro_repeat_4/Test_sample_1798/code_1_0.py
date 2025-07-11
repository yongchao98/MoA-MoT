import numpy as np
from scipy.integrate import quad

def solve():
    """
    Solves for alpha and P(X < 3) based on the problem description.
    """
    # Step 1 & 2: Define the circle parameters and the function f(x)
    # The equation of the circle passing through A(1,0), B(10,0), and C(4,3)
    # is (x - h)^2 + (y - k)^2 = r^2.
    h = 5.5
    k = -1.5
    r_squared = 22.5
    
    # The function y = f(x) is the upper semi-circle.
    # f(x) = k + sqrt(r_squared - (x-h)^2)
    def f(x):
        return k + np.sqrt(r_squared - (x - h)**2)

    # Step 3: Calculate alpha
    # The total area under f(x) from 1 to 10
    total_area, _ = quad(f, 1, 10)
    
    # alpha is the normalization constant
    alpha = 1 / total_area

    # Step 4: Calculate P(X < 3)
    # The area under f(x) from 1 to 3
    partial_area, _ = quad(f, 1, 3)
    
    # P(X < 3) is alpha * partial_area
    prob_X_lt_3 = alpha * partial_area
    
    print("The equation for the upper semi-circle is:")
    print(f"f(x) = {k} + sqrt({r_squared} - (x - {h})^2)")
    print("\nThe final values are calculated as follows:")
    print(f"alpha = 1 / (Integral of f(x) from 1 to 10) = 1 / {total_area} = {alpha}")
    print(f"P(X < 3) = alpha * (Integral of f(x) from 1 to 3) = {alpha} * {partial_area} = {prob_X_lt_3}")

solve()