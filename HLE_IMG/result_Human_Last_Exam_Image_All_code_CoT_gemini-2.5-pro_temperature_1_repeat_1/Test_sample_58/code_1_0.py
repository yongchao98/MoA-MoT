import sympy

# Define the symbolic variables used in the formula
h, r, theta, pi = sympy.symbols('h r theta pi')

# The problem asks for the height 'h' in terms of the radius 'r' and angle 'theta'.
# The derivation involves unrolling the cylinder and using geometric properties of the resulting paths.
# The final formula is constructed from several parts. We will print it piece by piece
# to show where the numbers in the equation come from, as requested.

print("The height h of the cylinder is given by the following equation:")
print("Note: The angle theta is in radians.")

# The formula can be built up as follows:
# h = r * ( (pi / 2) - theta ) * ( (1 + sin(theta)) / cos(theta) )
# The numbers in this equation are 2 and 1.

# Let's represent each component to construct the output string.
radius_term = 'r'
number_pi = 'pi'
number_2 = 2
angle_theta = 'theta'
number_1 = 1
sin_term = 'sin(theta)'
cos_term = 'cos(theta)'

# Print the equation in a way that highlights the numbers
print("\nHere is the equation with its components and numbers shown:")
print(f"h = {radius_term} * (({number_pi} / {number_2}) - {angle_theta}) * (({number_1} + {sin_term}) / {cos_term})")

# For a clearer mathematical representation, we use sympy's pretty printing feature.
print("\nIn standard mathematical notation, the formula is:")

# Define the expression for h using sympy
final_expression = r * (pi/2 - theta) * (1 + sympy.sin(theta)) / sympy.cos(theta)

# Create and print the final equation
final_equation = sympy.Eq(h, final_expression)
sympy.pprint(final_equation, use_unicode=True)
