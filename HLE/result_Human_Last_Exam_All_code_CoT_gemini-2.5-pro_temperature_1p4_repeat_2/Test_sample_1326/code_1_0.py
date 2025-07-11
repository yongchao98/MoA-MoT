import numpy as np
from scipy.integrate import quad

# The problem is to compute the integral of a piecewise function p(x) from 0 to 4.
# We must split the integral into two parts corresponding to the function's definition.

# Part 1: p(x) = (2*x^3)/8 for 0 <= x <= 3
def p_part1(x):
  """Defines the first piece of the function."""
  return (2 * x**3) / 8

# Part 2: p(x) = e^x * (1 + sin(x)) / (1 + cos(x)) for 3 < x <= 5
def p_part2(x):
  """Defines the second piece of the function."""
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# --- Calculations ---

# Calculate the integral for the first part from x=0 to x=3
integral_part1, error1 = quad(p_part1, 0, 3)

# Calculate the integral for the second part from x=3 to x=4
integral_part2, error2 = quad(p_part2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral_part1 + integral_part2

# --- Output ---

# Print the final equation with the calculated numerical values.
print("The total integral is the sum of the integrals of the two function pieces.")
print(f"Integral from 0 to 3: {integral_part1}")
print(f"Integral from 3 to 4: {integral_part2}")
print(f"Final Equation: {integral_part1} + ({integral_part2}) = {total_integral}")
