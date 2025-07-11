import math

# The problem is solved analytically to find an expression for r_0.
# The expression derived is r_0 = (226 + 49 * sqrt(2)) / 17.
# This script calculates the numerical value of this expression.

# The numbers from the final equation
a = 226
b = 49
c = 17
d = 2 # The number under the square root

# Calculate the value of r_0
r0 = (a + b * math.sqrt(d)) / c

# Print the equation and its calculated value
print("The radial distance r_0 is given by the equation:")
print(f"r_0 = ({a} + {b} * math.sqrt({d})) / {c}")
print("\nThe calculated value is:")
print(f"r_0 = {r0}")