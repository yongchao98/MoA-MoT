import sympy

# Define the symbolic variables for the equation
y, k, H, z1, beta = sympy.symbols('y k H z_1 beta')
psi = sympy.Function('psi')

# From the derivation, we have the slip-dependent zeta potential at the bottom wall
# Let's denote z_a1 = z1 * (1 + beta * k)
z_a1 = z1 * (1 + beta * k)

# Construct the right-hand side of the final derived solution based on the derivation:
# psi(y) = z_a1 * sinh(k*(H/2 - y)) / sinh(k*H)
numerator = z_a1 * sympy.sinh(k * (H/2 - y))
denominator = sympy.sinh(k * H)
rhs_expression = numerator / denominator

# Create the final equation for psi(y)
final_equation = sympy.Eq(psi(y), rhs_expression)

# Print the final expression
# The output will show each symbol from the final derived equation.
print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
print(final_equation)