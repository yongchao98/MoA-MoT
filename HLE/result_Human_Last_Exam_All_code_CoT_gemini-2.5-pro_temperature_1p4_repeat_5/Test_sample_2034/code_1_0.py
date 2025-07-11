import sympy

# Define symbols
c1, c2, beta, alpha_i = sympy.symbols('c1 c2 beta alpha_i')
k_alpha_i = 1 + beta * alpha_i

# The expression from the bound
bound_expr = (1 + c1*beta)*alpha_i - (1 + c2*beta)*k_alpha_i

# Based on the derivation, we determined the constants
c1_val = 2
c2_val = 0

final_bound_expr = bound_expr.subs({c1: c1_val, c2: c2_val})

# We are asked to return the final equation.
# The inequality is: - (K alpha_D-i)_i <= ... + o(beta)
# Let's print the determined RHS of the inequality.
# Let's use sympy to format it nicely.
final_rhs = sympy.simplify(final_bound_expr)
final_rhs_expanded = sympy.expand(final_rhs)

# We want to print the equation with all the numbers
# -(K alpha)_i <= (1 + 2*beta)*alpha_i - (1 + 0*beta)*(K alpha)_i
# -(K alpha)_i <= (1 + 2*beta)*alpha_i - (K alpha)_i

print("The extended Jaakola-Haussler bound is:")
print(f"- (K vec(alpha)^(D-i))_i <= (1 + {c1_val}*beta)*alpha_i - (1 + {c2_val}*beta)*(K vec(alpha)^D)_i + o(beta)")

# This simplifies to:
print("Simplified form:")
print(f"- (K vec(alpha)^(D-i))_i <= (1 + {c1_val}*beta)*alpha_i - (K vec(alpha)^D)_i + o(beta)")

# Let's present the determined values for c1 and c2.
print(f"\nDetermined values:")
print(f"c1 = {c1_val}")
print(f"c2 = {c2_val}")
