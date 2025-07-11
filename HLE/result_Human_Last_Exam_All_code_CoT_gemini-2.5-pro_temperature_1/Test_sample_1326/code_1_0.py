import numpy as np
from scipy.integrate import quad

# This script calculates the definite integral of a piecewise function from x=0 to x=4.

# Define the first part of the function for the interval [0, 3]
def p_part1(x):
  """
  Defines the function p(x) for 0 <= x <= 3.
  """
  return (2 * x**3) / 8

# Define the second part of the function for the interval [3, 5]
def p_part2(x):
  """
  Defines the function p(x) for 3 <= x <= 5.
  """
  # Note: np.cos(x) can be -1 at x = pi, which is approx 3.14159.
  # This value is within the integration range [3, 4], so we must handle potential division by zero.
  # However, quad can often handle singularities at the boundaries.
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Set the integration limits for the two parts
lower_bound1, mid_point, upper_bound2 = 0, 3, 4

# Calculate the integral for the first part
# The result of quad is a tuple (integral_value, estimated_error)
integral1, error1 = quad(p_part1, lower_bound1, mid_point)

# Calculate the integral for the second part
integral2, error2 = quad(p_part2, mid_point, upper_bound2)

# The total integral is the sum of the two parts
total_integral = integral1 + integral2

# Print the final result as an equation
print(f"The integral is split into two parts at x=3:")
print(f"Integral from {lower_bound1} to {mid_point}: {integral1}")
print(f"Integral from {mid_point} to {upper_bound2}: {integral2}")
print(f"Total Integral = {integral1} + ({integral2}) = {total_integral}")