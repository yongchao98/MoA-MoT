import sympy

# Define symbols for clarity
r_i, r_o = sympy.symbols('r_i r_o')

# The ratio of maximum tangential stress to internal pressure
ratio_formula = (r_o**2 + r_i**2) / (r_o**2 - r_i**2)

# The problem states the outer radius is twice the inner radius
# We can substitute r_o = 2*r_i into the formula.
# For calculation purposes, let's assume r_i = 1, which makes r_o = 2.
inner_radius_val = 1
outer_radius_val = 2

# Calculate the numerator and denominator of the ratio
numerator = outer_radius_val**2 + inner_radius_val**2
denominator = outer_radius_val**2 - inner_radius_val**2

# Calculate the final ratio
final_ratio = sympy.Rational(numerator, denominator)

# Print the explanation and the result
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
print("\nGiven the condition r_o = 2 * r_i, we can substitute values.")
print("Let r_i = 1, then r_o = 2.")
print("\nSubstituting these values into the formula:")
print(f"Ratio = ({outer_radius_val}^2 + {inner_radius_val}^2) / ({outer_radius_val}^2 - {inner_radius_val}^2)")
print(f"Ratio = ({outer_radius_val**2} + {inner_radius_val**2}) / ({outer_radius_val**2} - {inner_radius_val**2})")
print(f"Ratio = {numerator} / {denominator}")
print(f"\nThe final ratio is {final_ratio}.")
print(f"As a decimal, this is approximately {float(final_ratio):.3f}.")