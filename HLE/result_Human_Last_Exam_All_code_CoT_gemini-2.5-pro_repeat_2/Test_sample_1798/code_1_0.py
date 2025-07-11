import numpy as np
from scipy.integrate import quad

def solve_and_print_results():
    """
    Solves the problem by finding the circle equation, the normalization constant alpha,
    and the probability P(X < 3).
    """
    
    # Step 1 & 2: Define the circle parameters and the function f(x)
    # The circle passes through A=(1,0), B=(10,0), C=(4,3).
    # The center's x-coordinate (h) is the midpoint of AB's x-coordinates.
    h = (1 + 10) / 2
    
    # We use the general circle equation (x-h)^2 + (y-k)^2 = r^2
    # and substitute points A and C to find k and r^2.
    # For A(1,0): (1 - 5.5)^2 + (0 - k)^2 = r^2  => 20.25 + k^2 = r^2
    # For C(4,3): (4 - 5.5)^2 + (3 - k)^2 = r^2  => 2.25 + 9 - 6k + k^2 = r^2
    # Equating them: 20.25 = 11.25 - 6k => 6k = -9 => k = -1.5
    k = -1.5
    
    # Now find r^2
    r_squared = (1 - h)**2 + (0 - k)**2
    
    # The function f(x) is the upper semi-circle
    def f(x):
        return k + np.sqrt(r_squared - (x - h)**2)

    # Step 3: Find the normalization constant alpha
    # alpha = 1 / integral of f(x) from 1 to 10
    total_integral, _ = quad(f, 1, 10)
    alpha = 1 / total_integral

    # Step 4: Calculate the probability P(X < 3)
    # P(X < 3) = alpha * integral of f(x) from 1 to 3
    integral_lt_3, _ = quad(f, 1, 3)
    prob_lt_3 = alpha * integral_lt_3

    # Print the results as requested
    print("The equation of the circle is (x - h)^2 + (y - k)^2 = r^2, where:")
    print(f"h = {h}")
    print(f"k = {k}")
    print(f"r^2 = {r_squared}")
    print("-" * 20)
    print(f"The value of alpha is: {alpha}")
    print(f"The value of P(X < 3) is: {prob_lt_3}")

solve_and_print_results()