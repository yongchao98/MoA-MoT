import math

# The problem simplifies to calculating the value of the expression:
# (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
# This is derived by simplifying the denominator through a matrix optimization,
# and then decomposing the resulting integral into a sum of four parts,
# which are evaluated using standard integral forms related to the Gamma and Zeta functions.

# Constants from the derivation
pi = math.pi

# The coefficients and constant term in the final expression
coeff_pi8 = 8/15
coeff_pi2 = 1/3
coeff_pi = -1/2
const_term = 1

# Calculate the value of each part of the expression
part1 = coeff_pi8 * pi**8
part2 = coeff_pi2 * pi**2
part3 = coeff_pi * pi
part4 = const_term

# Calculate the final total value
total_value = part1 + part2 + part3 + part4

# As requested, we print each number in the final equation.
print("The integral evaluates to the following expression:")
print(f"({8}/{15}) * pi^8 + ({1}/{3}) * pi^2 - ({1}/{2}) * pi + {1}")
print("\nBreaking this down into its component values:")
print(f"Value of '(8/15) * pi^8': {part1}")
print(f"Value of '(1/3) * pi^2': {part2}")
print(f"Value of '(-1/2) * pi': {part3}")
print(f"Value of the constant term '1': {part4}")

print("\nThe final numerical result is:")
print(total_value)