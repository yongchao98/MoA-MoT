import numpy as np
from scipy.integrate import quad

# The function p(x) is defined in two parts.
# We define a function for each part.

# Part 1: p(x) = 2x^3 / 8 for 0 <= x <= 3
def p1(x):
  """
  This function represents the first part of p(x).
  """
  return (2 * x**3) / 8

# Part 2: p(x) = e^x(1 + sin(x)) / (1 + cos(x)) for 3 <= x <= 5
def p2(x):
  """
  This function represents the second part of p(x).
  """
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# The integral from 0 to 4 is split at x = 3.
# We calculate the integral for the first part from 0 to 3.
integral_part1, error1 = quad(p1, 0, 3)

# We calculate the integral for the second part from 3 to 4.
integral_part2, error2 = quad(p2, 3, 4)

# The total integral is the sum of the two parts.
total_integral = integral_part1 + integral_part2

# Print the results, showing each part of the calculation.
print("Integral of the first part (from 0 to 3):")
print(f"∫(2x^3/8)dx from 0 to 3 = {integral_part1}")
print("\nIntegral of the second part (from 3 to 4):")
print(f"∫(e^x(1+sinx)/(1+cosx))dx from 3 to 4 = {integral_part2}")
print("\nTotal integral from 0 to 4:")
print(f"Total = {integral_part1} + {integral_part2}")
print(f"Total = {total_integral}")
