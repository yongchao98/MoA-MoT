import numpy as np
from scipy import integrate

# Define the two parts of the piecewise function
def part1(x):
  """Defines the function for the interval 0 <= x <= 3."""
  return (2 * x**3) / 8

def part2(x):
  """Defines the function for the interval 3 <= x <= 5."""
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Calculate the integral of the first part from 0 to 3
# The result is a tuple (value, estimated_error)
integral_part1, error1 = integrate.quad(part1, 0, 3)

# Calculate the integral of the second part from 3 to 4
integral_part2, error2 = integrate.quad(part2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral_part1 + integral_part2

# Print the components of the final equation and the result
print("The integral is split into two parts: Integral_1 + Integral_2")
print(f"Integral_1 (from 0 to 3) = {integral_part1}")
print(f"Integral_2 (from 3 to 4) = {integral_part2}")
print(f"Total Integral = {integral_part1} + {integral_part2} = {total_integral}")
