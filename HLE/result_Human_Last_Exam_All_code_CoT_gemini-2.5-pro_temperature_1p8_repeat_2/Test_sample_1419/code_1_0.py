import math

# This script calculates and displays the leading-order expression for the
# Wilson-Fisher fixed point coupling (u*) in phi^4 theory near four dimensions.

# The derivation starts from the one-loop beta function in d = 4 - epsilon dimensions:
# β(u) = -ϵ*u + b*u^2

# The coefficient 'b' from the one-loop diagrammatic calculation is 3 / (16 * π^2).
coefficient_b_numerator = 3
coefficient_b_denominator = 16

# A fixed point u* is found by setting β(u*) = 0 and solving for the non-trivial solution:
# u* = ϵ / b

# Substituting the value of b gives:
# u* = ϵ / (3 / (16 * π^2))
# which simplifies to:
# u* = (16 * π^2 / 3) * ϵ

# Let's define the components of the final expression.
final_coefficient_numerator = 16
final_coefficient_denominator = 3
pi_symbol = "π"
epsilon_symbol = "ϵ"

print("The leading order expression for the fixed point coupling u* in the ϕ^4 theory near d=4 dimensions is:")
# We use an f-string to construct and print the final equation,
# ensuring each number is explicitly shown as requested.
print(f"u* = ({final_coefficient_numerator} * {pi_symbol}² / {final_coefficient_denominator}) * {epsilon_symbol}")

# For additional context, we calculate the approximate numerical value of the coefficient.
numerical_coefficient = (final_coefficient_numerator * math.pi**2) / final_coefficient_denominator
print(f"\nNumerically, the expression is approximately:")
print(f"u* ≈ {numerical_coefficient:.4f} * {epsilon_symbol}")