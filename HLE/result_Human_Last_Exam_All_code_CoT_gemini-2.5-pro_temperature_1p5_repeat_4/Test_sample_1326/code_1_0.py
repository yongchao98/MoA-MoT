import numpy as np
from scipy.integrate import quad

def solve_integral():
    """
    This function calculates the definite integral of a piecewise function p(x)
    from x = 0 to x = 4.

    The function p(x) is defined as:
    - p(x) = (2 * x**3) / 8                  if 0 <= x <= 3
    - p(x) = (e^x * (1 + sin x)) / (1 + cos x) if 3 < x <= 5
    """

    # Define the two parts of the function p(x)
    def p1(x):
        """p(x) for 0 <= x <= 3"""
        return (2 * x**3) / 8

    def p2(x):
        """p(x) for 3 < x <= 5"""
        # Using numpy for mathematical functions like exp, sin, and cos
        return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

    # Define the integration limits
    lower_bound = 0
    mid_point = 3
    upper_bound = 4

    # Calculate the integral for the first part (from 0 to 3)
    # The quad function returns the integral result and an error estimate.
    # We only need the result for this problem.
    integral_part1, _ = quad(p1, lower_bound, mid_point)

    # Calculate the integral for the second part (from 3 to 4)
    integral_part2, _ = quad(p2, mid_point, upper_bound)

    # Calculate the total integral by summing the two parts
    total_integral = integral_part1 + integral_part2

    # Print the final equation with the calculated values
    # The output format shows each number in the final equation.
    print(f"The integral from 0 to 4 is split into two parts at x=3.")
    print(f"The value for the integral from 0 to 3 is: {integral_part1:.4f}")
    print(f"The value for the integral from 3 to 4 is: {integral_part2:.4f}")
    print("\nThe final equation is:")
    print(f"{integral_part1:.4f} + ({integral_part2:.4f}) = {total_integral:.4f}")
    
    # Returning the final answer for automated checking
    return total_integral

# Execute the function to get the result
final_answer = solve_integral()
# The final answer is also printed here in the required format
# print(f"\n<<<{final_answer:.4f}>>>")