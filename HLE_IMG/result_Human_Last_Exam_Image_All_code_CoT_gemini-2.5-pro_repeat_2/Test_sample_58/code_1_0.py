import sympy

# Define the symbols for the variables
h, r, theta = sympy.symbols('h r theta')

# The relationship is derived from the geometry of the problem.
# We identified a right-angled triangle in 3D space with legs 'h' and '2*r'.
# The angle theta is related to this triangle.
# tan(theta) = opposite / adjacent = h / (2*r)
# From this, we solve for h.
equation = sympy.Eq(h, 2 * r * sympy.tan(theta))

# Print the equation for the height h
print("The height of the cylinder is given by the equation:")
# Using sympy.pretty_print for a nicer output format
sympy.pretty_print(equation, use_unicode=True)

# The prompt asks to output each number in the final equation.
# In the equation h = 2 * r * tan(theta), the only number is 2.
print("\nThe number in the final equation is: 2")
