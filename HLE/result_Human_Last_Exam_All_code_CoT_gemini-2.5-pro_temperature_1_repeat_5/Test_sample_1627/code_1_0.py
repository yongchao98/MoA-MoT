import fractions

# Problem statement: Find the ratio of maximum tangential stress to internal pressure
# for a thick-walled cylinder where the outer radius is twice the inner radius.

# Let the inner radius 'a' be 1 for simplicity.
a = 1
# The outer radius 'b' is twice the inner radius.
b = 2 * a

# The formula for the ratio is: (b^2 + a^2) / (b^2 - a^2)
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("\nGiven b = 2a, we can substitute these values into the formula.")
print(f"Let a = {a}, then b = {b}.")

# Calculate the terms
a_squared = a**2
b_squared = b**2

# Calculate the numerator and denominator of the ratio
numerator = b_squared + a_squared
denominator = b_squared - a_squared

# Calculate the final ratio
final_ratio = numerator / denominator

# Print the step-by-step calculation
print("\nSubstituting the values into the formula:")
print(f"Ratio = (({b})^2 + ({a})^2) / (({b})^2 - ({a})^2)")
print(f"Ratio = ({b_squared} + {a_squared}) / ({b_squared} - {a_squared})")
print(f"Ratio = {numerator} / {denominator}")

# Express the result as a fraction and a decimal
ratio_as_fraction = fractions.Fraction(numerator, denominator)

print(f"\nThe calculated ratio is {final_ratio:.4f}, which is equivalent to the fraction {ratio_as_fraction}.")