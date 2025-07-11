import numpy as np
from scipy.integrate import quad

def part1(x):
  """
  Defines the function p(x) for the interval 0 <= x <= 3.
  """
  return (2 * x**3) / 8

def part2(x):
  """
  Defines the function p(x) for the interval 3 <= x <= 5.
  """
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# The total integral must be split at x=3
# Calculate the integral of the first part from x=0 to x=3
integral_1, error_1 = quad(part1, 0, 3)

# Calculate the integral of the second part from x=3 to x=4
integral_2, error_2 = quad(part2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral_1 + integral_2

# Print the final result in an equation format
print(f"The integral from 0 to 4 is the sum of the integral from 0 to 3 and the integral from 3 to 4.")
print(f"Final Equation:")
print(f"{integral_1} + ({integral_2}) = {total_integral}")

# The following line is for the final answer extraction.
# The numeric value is the content of the <<<...>>> block.
# print(f"<<<{total_integral}>>>")