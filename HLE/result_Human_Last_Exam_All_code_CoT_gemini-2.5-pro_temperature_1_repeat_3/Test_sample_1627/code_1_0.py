import sympy

# Define the symbols for the radii
a, b = sympy.symbols('a b')

# The problem states that the outer radius is twice the inner radius.
# We can set a=1 and b=2 for the calculation, as the ratio is dimensionless.
inner_radius_val = 1
outer_radius_val = 2 * inner_radius_val

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)
numerator = outer_radius_val**2 + inner_radius_val**2
denominator = outer_radius_val**2 - inner_radius_val**2

# Calculate the final ratio
ratio = sympy.Rational(numerator, denominator)

# Print the step-by-step calculation
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("\nGiven the outer radius 'b' is twice the inner radius 'a', we can set a=1 and b=2.")
print("\nSubstituting these values into the formula:")
print(f"Ratio = ({outer_radius_val}^2 + {inner_radius_val}^2) / ({outer_radius_val}^2 - {inner_radius_val}^2)")
print(f"Ratio = ({outer_radius_val**2} + {inner_radius_val**2}) / ({outer_radius_val**2} - {inner_radius_val**2})")
print(f"Ratio = {numerator} / {denominator}")
print(f"\nThe final ratio is {ratio}.")
print(f"As a decimal, the ratio is approximately {float(ratio):.4f}.")