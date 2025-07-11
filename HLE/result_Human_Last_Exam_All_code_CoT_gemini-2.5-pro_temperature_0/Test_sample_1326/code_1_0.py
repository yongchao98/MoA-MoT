import numpy as np
from scipy.integrate import quad

# Define the first part of the function for the interval [0, 3]
def p1(x):
    """
    p(x) = 2*x^3 / 8 for 0 <= x <= 3
    """
    return (2 * x**3) / 8

# Define the second part of the function for the interval [3, 5]
def p2(x):
    """
    p(x) = e^x * (1 + sin(x)) / (1 + cos(x)) for 3 <= x <= 5
    """
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Calculate the integral of the first part from 0 to 3
# The result from quad is a tuple (value, error_estimate)
integral1, error1 = quad(p1, 0, 3)

# Calculate the integral of the second part from 3 to 4
integral2, error2 = quad(p2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral1 + integral2

# Print the results and the final equation
print("The integral is calculated by splitting it into two parts based on the function's definition.")
print(f"The integral from x=0 to x=3 is for the function p(x) = 2*x^3 / 8.")
print(f"The integral from x=3 to x=4 is for the function p(x) = e^x(1+sin(x))/(1+cos(x)).")
print("\n--- Calculation ---")
print(f"Value of the first integral (0 to 3): {integral1}")
print(f"Value of the second integral (3 to 4): {integral2}")
print("\nThe total integral is the sum of these two values.")
print(f"Final Equation: {integral1} + ({integral2}) = {total_integral}")
