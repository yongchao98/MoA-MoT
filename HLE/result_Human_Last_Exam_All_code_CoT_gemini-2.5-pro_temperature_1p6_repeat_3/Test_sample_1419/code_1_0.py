# This script calculates and displays the leading order expression for the
# Wilson-Fisher fixed point coupling in phi^4 theory near four dimensions.

# The derivation starts from the one-loop beta function in d = 4 - epsilon dimensions:
# beta(u) = -epsilon * u + (3 / (16 * pi^2)) * u^2
# The non-trivial fixed point u* is found by setting beta(u*) = 0.
# Solving for u* gives: u* = (16 * pi^2 / 3) * epsilon

# Define the numerical and symbolic components of the final expression.
numerator_coefficient = 16
denominator_coefficient = 3
pi_symbol = "π^2"
epsilon_symbol = "ε"

# Construct and print the final equation for the fixed point coupling.
# The result shows the leading order dependence on epsilon.
print("The leading order expression for the fixed point coupling u* in terms of ε is:")
print(f"u* = ( {numerator_coefficient} * {pi_symbol} / {denominator_coefficient} ) * {epsilon_symbol}")