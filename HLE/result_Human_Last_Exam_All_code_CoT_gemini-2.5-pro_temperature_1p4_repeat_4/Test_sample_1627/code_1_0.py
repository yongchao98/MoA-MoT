import fractions

# Define the relationship between outer radius (b) and inner radius (a)
# We are given that b = 2a.
# The actual value of 'a' doesn't matter for the ratio, so we can set a = 1 for simplicity.
a = 1
b = 2 * a

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)

# Calculate the terms in the equation
a_squared = a**2
b_squared = b**2

numerator = b_squared + a_squared
denominator = b_squared - a_squared

# Calculate the final ratio
ratio_decimal = numerator / denominator
ratio_fraction = fractions.Fraction(numerator, denominator)


# Print the steps and the final calculation
print("The formula for the ratio of maximum tangential stress to internal pressure is (b² + a²) / (b² - a²).")
print(f"Given b = 2a, we can set a = {a} and b = {b}.")
print("Substituting these values into the formula:")
print(f"Ratio = ({b}² + {a}²) / ({b}² - {a}²)")
print(f"Ratio = ({b_squared} + {a_squared}) / ({b_squared} - {denominator})")
print(f"Ratio = {numerator} / {denominator}")
print(f"The value of the ratio is {ratio_decimal}, which is equal to the fraction {ratio_fraction}.")