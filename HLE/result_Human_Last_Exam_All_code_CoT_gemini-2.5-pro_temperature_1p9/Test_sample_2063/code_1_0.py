# The final expression for the limiting average probability P(n) is (7/8)^8.
# This script calculates the value of this expression.

# The base of the fraction
numerator = 7
denominator = 8

# The exponent in the expression
exponent = 8

# Calculate the final probability
limiting_probability = (numerator / denominator) ** exponent

# Output the equation and the result, as requested.
print("Based on the derivation, the limiting probability is given by the expression:")
print(f"P(n -> infinity) = ({numerator}/{denominator})^{exponent}")
print("\nCalculating the components:")
print(f"Numerator: {numerator}^{exponent} = {numerator**exponent}")
print(f"Denominator: {denominator}^{exponent} = {denominator**exponent}")
print(f"\nFinal Limiting Value: {limiting_probability}")
