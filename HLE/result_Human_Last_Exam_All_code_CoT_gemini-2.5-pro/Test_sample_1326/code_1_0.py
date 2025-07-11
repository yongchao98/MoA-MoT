import numpy as np
from scipy.integrate import quad

# Define the first part of the function for the interval [0, 3]
def p1(x):
    """
    p(x) = 2x^3 / 8 for 0 <= x <= 3
    """
    return (2 * x**3) / 8

# Define the second part of the function for the interval [3, 5]
def p2(x):
    """
    p(x) = e^x(1 + sin x) / (1 + cos x) for 3 <= x <= 5
    """
    # Use np functions for vectorization and numerical stability
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Calculate the integral for the first part, from x = 0 to x = 3
integral1, error1 = quad(p1, 0, 3)

# Calculate the integral for the second part, from x = 3 to x = 4
integral2, error2 = quad(p2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral1 + integral2

# Print the final equation with the calculated values
print("The total integral is the sum of the integrals over the two intervals.")
print(f"Integral from 0 to 3 = {integral1}")
print(f"Integral from 3 to 4 = {integral2}")
print(f"Total Integral = {integral1} + ({integral2}) = {total_integral}")
