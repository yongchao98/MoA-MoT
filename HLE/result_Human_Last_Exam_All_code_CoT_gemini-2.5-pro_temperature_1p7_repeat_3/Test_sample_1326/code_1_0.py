import numpy as np
from scipy.integrate import quad

def p(x):
  """
  Defines the piecewise function p(x).
  """
  if 0 <= x <= 3:
    return (2 * x**3) / 8
  elif 3 < x <= 5:
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))
  else:
    return 0

# The integral needs to be split into two parts at x=3.

# Part 1: Integrate from 0 to 3
# For this interval, p(x) = (2x^3)/8
integral_part1, error1 = quad(lambda x: (2 * x**3) / 8, 0, 3)

# Part 2: Integrate from 3 to 4
# For this interval, p(x) = e^x * (1 + sin(x)) / (1 + cos(x))
integral_part2, error2 = quad(lambda x: (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x)), 3, 4)

# The total integral is the sum of the two parts
total_integral = integral_part1 + integral_part2

# Print the final equation with the calculated values
print(f"The integral is calculated by splitting it into two parts:")
print(f"Integral from 0 to 4 = (Integral from 0 to 3) + (Integral from 3 to 4)")
print(f"Integral = {integral_part1} + ({integral_part2})")
print(f"Total Integral = {total_integral}")
