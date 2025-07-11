import numpy as np
from scipy.integrate import quad

def p1(x):
  """
  Defines the function p(x) for the interval 0 <= x <= 3.
  """
  return (2 * x**3) / 8

def p2(x):
  """
  Defines the function p(x) for the interval 3 <= x <= 5.
  """
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# The integral is split at x=3
# Part 1: Integral from 0 to 3
integral1, error1 = quad(p1, 0, 3)

# Part 2: Integral from 3 to 4
integral2, error2 = quad(p2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral1 + integral2

print(f"The integral from x = 0 to x = 3 is calculated as: {integral1}")
print(f"The integral from x = 3 to x = 4 is calculated as: {integral2}")
print(f"The total integral from x = 0 to x = 4 is the sum of these two parts.")
print(f"Total Integral = {integral1} + ({integral2}) = {total_integral}")
