import numpy as np
from scipy.integrate import quad

# Define the first part of the function for the interval [0, 3]
def part1(x):
    return (2 * x**3) / 8

# Define the second part of the function for the interval [3, 5]
def part2(x):
    # Using np.exp, np.sin, np.cos for robust numerical computation
    return (np.exp(x) * (1 + np.sin(x))) / (1 + np.cos(x))

# Define the integration limits
lower_bound = 0
mid_point = 3
upper_bound = 4

# Calculate the integral for the first part [0, 3]
# The result is a tuple (value, estimated_error)
integral1, error1 = quad(part1, lower_bound, mid_point)

# Calculate the integral for the second part [3, 4]
integral2, error2 = quad(part2, mid_point, upper_bound)

# The total integral is the sum of the two parts
total_integral = integral1 + integral2

# Print the value of each integral and the final equation
print(f"The integral from {lower_bound} to {mid_point} is: {integral1}")
print(f"The integral from {mid_point} to {upper_bound} is: {integral2}")
print("\nThe total integral is the sum:")
print(f"{integral1} + ({integral2}) = {total_integral}")