import fractions

# Problem Statement:
# Find the ratio of the maximum tangential stress to the internal pressure for a
# thick-walled cylinder where the outer radius is twice its inner radius.

# Let 'a' be the inner radius and 'b' be the outer radius.
# According to the problem, b = 2a.
# For simplicity in calculation, we can assume a = 1, which makes b = 2.
a = 1
b = 2

# The formula for the ratio of maximum tangential stress (at r=a) to the
# internal pressure (P_i) is given by:
# Ratio = (b^2 + a^2) / (b^2 - a^2)

# First, calculate the squared values of the radii.
a_squared = a**2
b_squared = b**2

# Now, calculate the numerator and the denominator of the formula.
numerator = b_squared + a_squared
denominator = b_squared - a_squared

# Finally, calculate the ratio.
ratio_value = numerator / denominator

# To present the result clearly, we show the equation with the numbers substituted.
# The prompt requires outputting each number in the final equation.
print(f"The formula for the ratio is: (b² + a²) / (b² - a²)")
print(f"Given b = 2a, we can set a = {a} and b = {b}.")
print("Substituting these values into the equation:")
print(f"Ratio = (({b}²) + ({a}²)) / (({b}²) - ({a}²))")
print(f"Ratio = ({b_squared} + {a_squared}) / ({b_squared} - {a_squared})")
print(f"Ratio = {numerator} / {denominator}")

# The exact ratio can be expressed as a fraction.
ratio_fraction = fractions.Fraction(numerator, denominator)
print(f"\nThe exact ratio as a fraction is: {ratio_fraction.numerator}/{ratio_fraction.denominator}")
# This corresponds to answer choice E.

# The ratio as a decimal is:
print(f"The ratio as a decimal is approximately: {ratio_value:.3f}")
# This corresponds to answer choice H.