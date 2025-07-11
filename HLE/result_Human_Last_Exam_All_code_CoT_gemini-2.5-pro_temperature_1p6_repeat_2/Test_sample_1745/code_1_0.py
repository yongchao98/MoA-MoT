import sympy

# Define the symbolic variables used in the equation
y, k, H, z1, beta = sympy.symbols('y k H z_1 beta')

# Construct the final expression for the potential distribution psi(y)
# This is derived from solving the linearized Poisson-Boltzmann equation
# with boundary conditions psi(H/2) = 0 and psi(-H/2) = z1*(1 + beta*k).
slip_dependent_zeta = z1 * (1 + beta * k)
numerator = sympy.sinh(k * (H/2 - y))
denominator = sympy.sinh(k * H)

psi_expression = slip_dependent_zeta * numerator / denominator

# Print the final resulting equation in a formatted way
print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
sympy.pprint(psi_expression, use_unicode=True)

# As requested, explicitly mentioning the numbers in the final equation.
print("\nBreakdown of the final equation to highlight the numbers:")
print("The full equation is: psi(y) = z_1*(1 + beta*k) * sinh(k*(H/2 - y)) / sinh(k*H)")
print("\nThe numbers present in this equation are:")
print("1: Found in the slip-dependent zeta potential term (1 + beta*k).")
print("1/2: Found in the numerator's argument as H/2, representing half the channel height.")
