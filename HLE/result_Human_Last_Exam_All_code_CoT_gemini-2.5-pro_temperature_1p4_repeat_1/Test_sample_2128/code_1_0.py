import math

# The problem asks for 1/p_n for n = 1000.
n = 1000

# The formula for 1/p_n is 4 * cos^2(pi / (n + 2)).
# Here, the numbers in the final equation are 4 and (n+2).
# Let's define them.
constant_multiplier = 4
denominator_in_angle = n + 2

# Calculate the final value
result = constant_multiplier * (math.cos(math.pi / denominator_in_angle))**2

# Print the final equation with the numbers and the result
print(f"The value is calculated using the formula: {constant_multiplier} * cos(pi / {denominator_in_angle})^2")
print(f"The final result is: {result}")