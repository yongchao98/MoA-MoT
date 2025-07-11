import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    Calculates the definite integral of a piecewise function p(x) from x=0 to x=4.
    """
    # Define the first part of the function for the interval [0, 3]
    def p1(x):
        return (2 * x**3) / 8

    # Define the second part of the function for the interval [3, 5]
    def p2(x):
        # Use np functions for compatibility with scipy.integrate.quad
        return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

    # Calculate the integral for the first part, from x = 0 to x = 3
    integral_part1, error1 = quad(p1, 0, 3)

    # Calculate the integral for the second part, from x = 3 to x = 4
    integral_part2, error2 = quad(p2, 3, 4)

    # The total integral is the sum of the two parts
    total_integral = integral_part1 + integral_part2

    # Print the final equation with the calculated values
    print("The integral is calculated by splitting it into two parts based on the function definition:")
    print(f"Total Integral = (Integral from 0 to 3) + (Integral from 3 to 4)")
    print(f"Total Integral = {integral_part1:.4f} + ({integral_part2:.4f})")
    print(f"Total Integral = {total_integral:.4f}")

solve_integral()