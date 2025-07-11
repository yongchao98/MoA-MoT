import math

# Define the constants from the problem
base = 12
power = 4

# Based on the analytical solution, the value of the integral expression is:
# I = (5^(1/4)) / 12
integral_value = (5**0.25) / 12

# Calculate the final result
final_result = (base**power) * (integral_value**power)

# Print the final equation with each number.
# The result is an integer, so we format it as such.
print(f"({base})^{power} * ({integral_value})^{power} = {final_result:.0f}")