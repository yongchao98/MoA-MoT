import numpy as np
from scipy.optimize import fsolve

# Define the transcendental equation to be solved
def equation_to_solve(x):
    """
    Represents the equation tan(x) - x - pi/2 = 0.
    """
    return np.tan(x) - x - np.pi / 2

# Find the root of the equation.
# From graphical analysis, the root is in the interval (0, pi/2).
# We use an initial guess of 1.3.
x0 = fsolve(equation_to_solve, 1.3)[0]

# Calculate the constant C using the solution x0
numerator = 2 * x0 + np.pi + 2 * np.sin(x0)
denominator = 2 * x0 + np.pi - 2 * np.sin(x0)
C = numerator / denominator

# Print the components of the final calculation
print(f"The optimal value x is the solution to tan(x) = x + pi/2.")
print(f"Numerically, x = {x0}")
print(f"The constant C is given by the expression (2*x + pi + 2*sin(x)) / (2*x + pi - 2*sin(x)).")
print(f"Plugging in the values, we get:")
print(f"Numerator = 2 * {x0} + {np.pi} + 2 * {np.sin(x0)} = {numerator}")
print(f"Denominator = 2 * {x0} + {np.pi} - 2 * {np.sin(x0)} = {denominator}")
print(f"C = {numerator} / {denominator}")
print(f"The smallest possible constant C is: {C}")
