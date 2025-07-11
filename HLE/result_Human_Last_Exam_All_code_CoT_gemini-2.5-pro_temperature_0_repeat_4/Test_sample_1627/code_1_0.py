import fractions

# Define the relationship between the inner radius 'a' and outer radius 'b'.
# The problem states that the outer radius is twice the inner radius (b = 2a).
# For calculation purposes, we can set a = 1, which means b = 2.
a = 1
b = 2

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)

# Calculate the terms needed for the formula
a_squared = a**2
b_squared = b**2
numerator = b_squared + a_squared
denominator = b_squared - a_squared

# Calculate the final ratio
ratio_decimal = numerator / denominator
ratio_fraction = fractions.Fraction(numerator, denominator)

# Print the explanation and the step-by-step calculation
print("The formula for the ratio of maximum tangential stress to internal pressure is (b^2 + a^2) / (b^2 - a^2).")
print(f"Given that the outer radius 'b' is twice the inner radius 'a', we can set a = {a} and b = {b}.")
print("\nSubstituting these values into the formula:")
print(f"Ratio = ({b}^2 + {a}^2) / ({b}^2 - {a}^2)")
print(f"Ratio = ({b_squared} + {a_squared}) / ({b_squared} - {a_squared})")
print(f"Ratio = {numerator} / {denominator}")
print(f"\nThe final ratio is {ratio_decimal}, which is equivalent to the fraction {ratio_fraction}.")
