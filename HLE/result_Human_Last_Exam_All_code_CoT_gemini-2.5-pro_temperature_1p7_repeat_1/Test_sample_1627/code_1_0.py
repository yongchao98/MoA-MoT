import sympy

# Problem Statement: Find the ratio of maximum tangential stress to internal pressure
# for a thick-walled cylinder where the outer radius (b) is twice the inner radius (a).

# The formula for the ratio is (b^2 + a^2) / (b^2 - a^2).
# We are given the condition b = 2a.

# Let's represent the ratio of radii
b_over_a = 2

# The formula can be simplified by dividing numerator and denominator by a^2:
# Ratio = ( (b/a)^2 + 1 ) / ( (b/a)^2 - 1 )

# Substitute the value of b/a = 2 into the simplified formula.
numerator_term_1 = b_over_a**2
denominator_term_1 = b_over_a**2

numerator = numerator_term_1 + 1
denominator = denominator_term_1 - 1

# Using sympy to display the fraction neatly
ratio_fraction = sympy.Rational(numerator, denominator)

print(f"The formula for the ratio of maximum tangential stress to internal pressure is (b^2 + a^2) / (b^2 - a^2).")
print(f"Given that the outer radius is twice the inner radius (b = 2a), the ratio b/a is {b_over_a}.")
print(f"Substituting this into the formula gives: ({b_over_a}^2 + 1) / ({b_over_a}^2 - 1)")
print(f"Calculating the squared terms: ({numerator_term_1} + 1) / ({denominator_term_1} - 1)")
print(f"The final equation is: {numerator} / {denominator}")
print(f"The result is {ratio_fraction}, which is approximately {float(ratio_fraction):.4f}.")