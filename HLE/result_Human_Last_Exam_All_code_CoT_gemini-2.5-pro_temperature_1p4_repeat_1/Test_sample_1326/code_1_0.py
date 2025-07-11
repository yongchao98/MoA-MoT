import numpy as np
from scipy.integrate import quad

# Define the first part of the function for the range [0, 3]
def p1(x):
  """
  Represents p(x) for 0 <= x <= 3
  """
  return (2 * x**3) / 8

# Define the second part of the function for the range [3, 5]
def p2(x):
  """
  Represents p(x) for 3 <= x <= 5
  """
  return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Calculate the integral of the first part from 0 to 3
# The quad function returns the integral and an estimated error
integral1, error1 = quad(p1, 0, 3)

# Calculate the integral of the second part from 3 to 4
integral2, error2 = quad(p2, 3, 4)

# The total integral is the sum of the two parts
total_integral = integral1 + integral2

# Print the results in a formatted equation
print(f"The integral from 0 to 3 is: {integral1}")
print(f"The integral from 3 to 4 is: {integral2}")
print(f"The total integral from 0 to 4 is: {integral1} + ({integral2}) = {total_integral}")