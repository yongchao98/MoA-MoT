import sympy

# Define the symbols for the radii
r_i, r_o = sympy.symbols('r_i r_o')

# The ratio of max tangential stress to internal pressure is given by the formula:
ratio_formula = (r_o**2 + r_i**2) / (r_o**2 - r_i**2)

# The problem states that the outer radius is twice the inner radius
# We can substitute r_o = 2*r_i into the formula.
# For calculation purposes, let's set r_i = 1, which means r_o = 2.
inner_radius = 1
outer_radius = 2

# Calculate the numerator and denominator of the fraction
numerator_val = outer_radius**2 + inner_radius**2
denominator_val = outer_radius**2 - inner_radius**2

# Calculate the final ratio
result = numerator_val / denominator_val

# Display the steps of the calculation
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
print("\nGiven r_o = 2 * r_i, we can substitute values like r_o = 2 and r_i = 1.")
print("\nSubstituting these values into the formula:")
print(f"Ratio = ({outer_radius}^2 + {inner_radius}^2) / ({outer_radius}^2 - {inner_radius}^2)")
print(f"Ratio = ({outer_radius**2} + {inner_radius**2}) / ({outer_radius**2} - {inner_radius**2})")
print(f"Ratio = {numerator_val} / {denominator_val}")
print(f"\nThe numerical value of the ratio is {result:.4f}, which as a fraction is {numerator_val}/{denominator_val}.")
