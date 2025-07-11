import sympy as sp

# Define the symbols
r, theta, h = sp.symbols('r theta h')

# The relationship derived from a standard interpretation of this problem geometry is
# h = 2 * r * tan(theta).
# Here, we will express h in terms of r and theta and print the relationship.

equation = sp.Eq(h, 2 * r * sp.tan(theta))

# Print the relationship
print("The height 'h' of the cylinder is given by the equation:")
print(f"{sp.pretty(equation.lhs)} = {sp.pretty(equation.rhs)}")
