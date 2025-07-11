import sympy

# Define the symbols for the radii
r_i, r_o = sympy.symbols('r_i r_o')

# The formula for the ratio of maximum tangential stress to internal pressure
ratio_formula = (r_o**2 + r_i**2) / (r_o**2 - r_i**2)

# The given condition is that the outer radius is twice the inner radius
# We can substitute r_o = 2*r_i into the formula
# To make it a numerical calculation, let's assume r_i = 1, which means r_o = 2.
inner_radius = 1
outer_radius = 2

# Calculate the numerator and denominator
numerator = outer_radius**2 + inner_radius**2
denominator = outer_radius**2 - inner_radius**2

# Calculate the final ratio
result = numerator / denominator

# Print the step-by-step calculation
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)\n")
print(f"Given that the outer radius is twice the inner radius (r_o = 2 * r_i).")
print(f"Let's set the inner radius r_i = {inner_radius} and the outer radius r_o = {outer_radius}.\n")
print("Substituting these values into the formula:")
print(f"Ratio = ({outer_radius}^2 + {inner_radius}^2) / ({outer_radius}^2 - {inner_radius}^2)")
print(f"Ratio = ({outer_radius**2} + {inner_radius**2}) / ({outer_radius**2} - {inner_radius**2})")
print(f"Ratio = {numerator} / {denominator}\n")
print(f"The calculated ratio is {numerator}/{denominator}, which is approximately {result:.4f}.")