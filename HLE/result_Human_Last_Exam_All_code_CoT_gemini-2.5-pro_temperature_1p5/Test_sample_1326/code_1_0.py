import numpy as np
from scipy.integrate import quad

# Define the first part of the function for the interval [0, 3]
def p_first_part(x):
    """
    Function definition for 0 <= x <= 3.
    p(x) = 2*x**3 / 8
    """
    return (x**3) / 4

# Define the second part of the function for the interval [3, 5]
def p_second_part(x):
    """
    Function definition for 3 <= x <= 5.
    p(x) = e**x * (1 + sin(x)) / (1 + cos(x))
    """
    # Note: np.sin and np.cos work with radians, which is standard for calculus.
    return np.exp(x) * (1 + np.sin(x)) / (1 + np.cos(x))

# Calculate the integral for the first part, from x = 0 to x = 3
integral_part1, error_part1 = quad(p_first_part, 0, 3)

# Calculate the integral for the second part, from x = 3 to x = 4
integral_part2, error_part2 = quad(p_second_part, 3, 4)

# Calculate the total integral
total_integral = integral_part1 + integral_part2

# Print the results in the required format
print(f"The integral is calculated by splitting it into two parts:")
print(f"Integral from 0 to 3 of (x^3 / 4) dx = {integral_part1}")
print(f"Integral from 3 to 4 of (e^x * (1 + sin(x)) / (1 + cos(x))) dx = {integral_part2}")
print(f"\nThe total integral from 0 to 4 is the sum of these two parts:")
print(f"Total Integral = {integral_part1} + ({integral_part2}) = {total_integral}")
