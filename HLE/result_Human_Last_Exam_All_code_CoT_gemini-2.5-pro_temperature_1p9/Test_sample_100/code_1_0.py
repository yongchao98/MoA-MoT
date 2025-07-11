import math

# Based on the detailed derivation, the integral evaluates to the expression:
# (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
# This script calculates the numerical value of this expression.
# The instruction "output each number in the final equation" is fulfilled by
# printing the full expression before displaying the final numerical result.

pi = math.pi

# The coefficients and constants from the final derived expression
c1 = 8/15
c2 = 1/3
c3 = -1/2
c4 = 1

# Calculate the final result
result = c1 * (pi**8) + c2 * (pi**2) + c3 * pi + c4

# Print the final equation with all its numeric components and the calculated value.
# The numbers are 8, 15, 1, 3, -1, 2, 1 which are all present in the expression below.
print(f"The integral evaluates to the expression: 8/15 * pi**8 + 1/3 * pi**2 - 1/2 * pi + 1")
print(f"Calculating the final value:")
print(f"8/15 * ({pi})**8 + 1/3 * ({pi})**2 - 1/2 * ({pi}) + 1 = {result}")