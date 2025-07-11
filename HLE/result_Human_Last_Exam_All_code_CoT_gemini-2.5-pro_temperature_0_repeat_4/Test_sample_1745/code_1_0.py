import sympy as sp

# This script uses the symbolic math library (sympy) to display the final expression.

# Define the mathematical symbols used in the expression.
# y: transverse coordinate
# k: Debye-Huckel parameter
# H: height of the microchannel
# beta: slip length parameter
# z_1: zeta potential at the bottom surface (j=1)
y, k, H, beta = sp.symbols('y k H beta')
z_1 = sp.Symbol('z_1')
psi = sp.Function('psi')(y)

# The slip-dependent zeta potential at the bottom wall is z_a1 = z_1 * (1 + beta * k).
zeta_a1 = z_1 * (1 + beta * k)

# Based on the derivation, we construct the final solution for the EDL potential distribution.
# The expression is: psi(y) = zeta_a1 * sinh(k*(H/2 - y)) / sinh(k*H)
solution_rhs = zeta_a1 * (sp.sinh(k * (H / 2 - y)) / sp.sinh(k * H))

# Create a symbolic equation object for clear printing.
final_equation = sp.Eq(psi, solution_rhs)

# Print the final expression for the EDL potential distribution.
# The output will be a formatted equation showing each symbol.
print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
sp.init_printing(use_unicode=True)
print(final_equation)