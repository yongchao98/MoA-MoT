import numpy as np
from scipy.integrate import quad

def solve_and_calculate():
    """
    This function solves the problem by finding the circle equation,
    calculating alpha, and then calculating the probability P(X < 3).
    """
    # Step 1 & 2: Determine the circle equation and the function f(x)
    # From solving the system of equations for points A(1,0), B(10,0), C(4,3):
    # (1-h)^2 + k^2 = (10-h)^2 + k^2  => 1 - 2h = 100 - 20h => 18h = 99 => h = 5.5
    # Using h=5.5 in the equations for A and C:
    # (1-5.5)^2 + k^2 = r^2 => 20.25 + k^2 = r^2
    # (4-5.5)^2 + (3-k)^2 = r^2 => 2.25 + 9 - 6k + k^2 = r^2
    # 20.25 + k^2 = 11.25 - 6k + k^2 => 9 = -6k => k = -1.5
    # r^2 = 20.25 + (-1.5)^2 = 20.25 + 2.25 = 22.5
    h = 5.5
    k = -1.5
    r_squared = 22.5

    # The function f(x) is the upper semi-circle equation
    def f(x):
        radicand = r_squared - (x - h)**2
        # The integration range [1, 10] is within the circle's domain on the x-axis,
        # so radicand will not be negative.
        return k + np.sqrt(radicand)

    print("The equation of the circle is (x - h)^2 + (y - k)^2 = r^2")
    print(f"h = {h}")
    print(f"k = {k}")
    print(f"r^2 = {r_squared}")
    print("-" * 20)

    # Step 3: Calculate alpha
    # alpha = 1 / integral of f(x) from 1 to 10
    total_integral, total_integral_error = quad(f, 1, 10)
    
    if total_integral <= 0:
        print("Error: The integral of f(x) must be positive.")
        return

    alpha = 1 / total_integral

    print(f"The value of alpha is: {alpha}")

    # Step 4: Calculate P(X < 3)
    # P(X < 3) = alpha * integral of f(x) from 1 to 3
    partial_integral, partial_integral_error = quad(f, 1, 3)
    
    prob_X_less_than_3 = alpha * partial_integral

    print(f"The value of P(X < 3) is: {prob_X_less_than_3}")

    # For the final answer format
    # print(f"<<<{alpha},{prob_X_less_than_3}>>>")


solve_and_calculate()