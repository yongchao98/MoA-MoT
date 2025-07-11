import math

# This script derives and prints the expression for the Wilson-Fisher fixed point
# in phi^4 theory near d=4 dimensions.

# The one-loop beta function is of the form:
# beta(u) = -epsilon * u + B * u^2
# For phi^4 theory, the coefficient B is 3 / (16 * pi^2).

# A fixed point u* satisfies beta(u*) = 0:
# -epsilon * u* + B * (u*)^2 = 0
# For the non-trivial fixed point (u* != 0), we can divide by u*:
# -epsilon + B * u* = 0
# This gives u* = epsilon / B

# Substitute the value of B:
# u* = epsilon / (3 / (16 * pi^2))
# u* = (16 * pi^2 / 3) * epsilon

# Define the components of the final expression for printing
coupling_symbol = "u*"
epsilon_symbol = "ε"
pi_symbol = "π"
numerator = 16
denominator = 3

# Print the final expression step-by-step
print(f"The beta function for the coupling u is: β(u) = -{epsilon_symbol}*u + ({denominator} / ({numerator}*{pi_symbol}²)) * u²")
print(f"Setting β({coupling_symbol}) = 0 gives the fixed point.")
print(f"The non-trivial solution is:")
print(f"{coupling_symbol} = ({numerator} * {pi_symbol}² / {denominator}) * {epsilon_symbol}")