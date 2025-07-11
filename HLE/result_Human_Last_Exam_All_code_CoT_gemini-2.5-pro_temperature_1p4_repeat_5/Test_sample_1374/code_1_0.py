import math

# This script calculates the furthest distance from a point A to the surface of a
# 1 m³ body of uniform density, shaped to produce the maximum possible gravitational
# field at A.

# From the physics of the problem, this distance is given by the formula:
# d = (15 / (4 * pi))^(1/3)

# Define the numerical constants from the formula.
numerator = 15
denominator_coeff = 4
pi_val = math.pi
exponent_val = 1/3

# Calculate the final result.
result = (numerator / (denominator_coeff * pi_val)) ** exponent_val

# Print the equation with each number explicitly shown, along with the final answer.
# This fulfills the request to "output each number in the final equation".
print(
    f"The furthest distance is calculated by the equation: "
    f"({numerator} / ({denominator_coeff} * {pi_val})) ^ ({exponent_val}) "
    f"= {result}"
)
print(f"\nThe furthest point on the surface of the playdough is {result:.5f} meters from point A.")
