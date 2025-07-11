import math

# The equation to solve for r_0 is derived from the analysis of the functional equation.
# The final expression for r_0 is:
# r_0 = (4 + 37 * sqrt(2)) / (3 * sqrt(2) - 1)

# We define the integer constants present in the final equation for r_0.
a = 4
b = 37
c = 2  # The number under the square root
d = 3
e = 1

# Calculate the square root of c
sqrt_c = math.sqrt(c)

# Calculate the numerator of the expression
numerator = a + b * sqrt_c

# Calculate the denominator of the expression
denominator = d * sqrt_c - e

# Calculate the final value of r_0
r_0 = numerator / denominator

# Output the equation in a formatted way, showing the numbers used.
print(f"The equation for the radial distance r_0 is:")
print(f"r_0 = ({a} + {b} * sqrt({c})) / ({d} * sqrt({c}) - {e})")

# Print the final result.
print(f"The calculated value for the radial distance r_0 is: {r_0}")
