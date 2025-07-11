import sympy

# This script generates and prints the symbolic expression for the
# Electrical Double-Layer (EDL) potential distribution, psi(y).

# 1. Define all the symbolic variables used in the equation.
y, H, k, beta, z1 = sympy.symbols('y H k beta z_1')
psi = sympy.Function('psi')(y)

# 2. Construct the right-hand side of the final derived equation.
# The number '1' is explicitly included in the term (1 + beta*k).
rhs = z1 * (1 + beta * k) * (sympy.sinh(k * (H - y)) / sympy.sinh(k * H))

# 3. Create the full equation object for pretty printing.
final_equation = sympy.Eq(psi, rhs)

# 4. Print the final equation in a readable format.
# The printed result includes all variables and numbers as they appear in the formula.
print("The final expression for the Electrical double-layer potential distribution psi(y) is:")
sympy.pprint(final_equation, use_unicode=True)