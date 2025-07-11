import numpy as np
from scipy import integrate

def p1(x):
  """
  Defines the first part of the function p(x) for 0 <= x <= 3.
  """
  return (2 * x**3) / 8

def p2(x):
  """
  Defines the second part of the function p(x) for 3 < x <= 5.
  """
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# The integration range is from 0 to 4. We split this at x = 3.
# The split point where the function definition changes.
split_point = 3
# The overall integration bounds.
lower_bound = 0
upper_bound = 4

# Calculate the integral for the first part: from lower_bound to split_point.
integral1, error1 = integrate.quad(p1, lower_bound, split_point)

# Calculate the integral for the second part: from split_point to upper_bound.
integral2, error2 = integrate.quad(p2, split_point, upper_bound)

# The total integral is the sum of the two parts.
total_integral = integral1 + integral2

# Print the results in a clear, step-by-step manner.
print(f"The integral is split at x = {split_point} due to the piecewise nature of the function.")
print(f"First integral (from x={lower_bound} to x={split_point}): {integral1}")
print(f"Second integral (from x={split_point} to x={upper_bound}): {integral2}")
print("\nThe final result is the sum of these two integrals.")
print("Final Equation:")
print(f"{integral1} + ({integral2}) = {total_integral}")

# Print the final numerical answer in the specified format
# print(f"\nFinal Answer: {total_integral}")