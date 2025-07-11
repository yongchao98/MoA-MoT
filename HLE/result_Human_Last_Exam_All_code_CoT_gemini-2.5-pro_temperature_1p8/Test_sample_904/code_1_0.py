import sympy

# Define the symbols used in the problem
# gamma is the surface tension
# r is the radial position
gamma, r = sympy.symbols('γ r')

# The governing equation is of the form:
# A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0

# From the derivation based on the linearized Young-Laplace equation,
# the coefficients A(r) and B(r) are determined by the surface tension term.
A_r = gamma
B_r = gamma / r

# The term C(r, ξ) corresponds to the external forces, which in this
# case is the electrostatic pressure, but its explicit form is not
# needed to find A(r) and B(r).

# Print the final expressions for A(r) and B(r)
print("Based on the derivation from the balance of pressures at the fluid interface:")
print("\nThe coefficient A(r) in the governing equation is:")
print(f"A(r) = {A_r}")

print("\nThe coefficient B(r) in the governing equation is:")
print(f"B(r) = {B_r}")

# To demonstrate the final equation structure:
d2xi_dr2 = sympy.Function('d²ξ/dr²')
dxi_dr = sympy.Function('dξ/dr')
C_term = sympy.Function('C')(r, sympy.Function('ξ')(r))

# Print the full equation structure
print("\nThe full linear equation for the interfacial shape ξ(r) is:")
equation = sympy.Eq(A_r * d2xi_dr2 + B_r * dxi_dr + C_term, 0)
# Using sympy.pprint for a more readable mathematical output
sympy.pprint(equation, use_unicode=True)