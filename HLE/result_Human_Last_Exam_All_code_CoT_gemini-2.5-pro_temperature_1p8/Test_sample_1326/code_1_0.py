import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Calculates the integral of the piecewise function p(x) from x = 0 to x = 4.
    """
    # Define the first part of the function for the interval [0, 3]
    def p1(x):
        return (2 * x**3) / 8

    # Define the second part of the function for the interval [3, 5]
    def p2(x):
        return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

    # 1. Calculate the integral for the first part: from x = 0 to x = 3
    integral1, _ = quad(p1, 0, 3)

    # 2. Calculate the integral for the second part: from x = 3 to x = 4
    # The function p2 has a singularity at x = pi because the denominator 1 + cos(x) becomes 0.
    # We pass this point to the 'quad' function for a more accurate result.
    singularity_point = np.pi
    integral2, _ = quad(p2, 3, 4, points=[singularity_point])

    # 3. Sum the two integrals to get the total integral
    total_integral = integral1 + integral2

    # 4. Print the final equation with the calculated values
    print("The total integral is the sum of the integrals of the two parts of the function.")
    print(f"First part (from 0 to 3): {integral1:.4f}")
    print(f"Second part (from 3 to 4): {integral2:.4f}")
    print("\nFinal Equation:")
    print(f"{integral1:.4f} + ({integral2:.4f}) = {total_integral:.4f}")

solve_integral()