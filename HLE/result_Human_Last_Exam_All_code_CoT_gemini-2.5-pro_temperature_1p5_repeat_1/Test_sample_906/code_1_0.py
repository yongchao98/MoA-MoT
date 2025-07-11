import math

# The problem asks for the symbolic expression for the steady-state probability pi_0.
# The derivation, as explained above, leads to pi_0 = e^(-rho).
# This code will print this final expression.
# The components of the equation are the variable 'pi_0', the equals sign,
# the base of the natural logarithm 'e', an exponentiation operator,
# a negative sign, and the variable 'rho'.

equation_lhs = "pi_0"
equals_sign = "="
expression_rhs = "e**(-rho)"

print(f"The final expression for the steady-state probability is:")
print(f"{equation_lhs} {equals_sign} {expression_rhs}")