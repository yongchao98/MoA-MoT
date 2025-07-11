import sympy as sp

# Define symbols
c1 = sp.Symbol('c1')
c2 = sp.Symbol('c2')
beta = sp.Symbol('beta')
alpha_i_D = sp.Symbol('alpha_i^D')
K_alpha_D_i = sp.Symbol('(K*alpha^D)_i')
K_alpha_D_minus_i_i = sp.Symbol('(K*alpha^(D-i))_i')
o_beta = sp.Symbol('o(beta)')

# Set the determined values
c1_val = -1
c2_val = 1

# Construct the inequality expression
# LHS: - (K * alpha^{D-i})_i
# RHS: +(1 + c1 * beta) * alpha_i^D - (1 + c2 * beta) * (K * alpha^D)_i + o(beta)
rhs_expression = (1 + c1 * beta) * alpha_i_D - (1 + c2 * beta) * K_alpha_D_i + o_beta
inequality = sp.Le(-K_alpha_D_minus_i_i, rhs_expression)

# Substitute the determined values of c1 and c2
final_inequality = inequality.subs({c1: c1_val, c2: c2_val})

# Print the determined values and the final equation
print(f"The determined values are:")
print(f"c1 = {c1_val}")
print(f"c2 = {c2_val}\n")

print("The extended Jaakola-Haussler bound is:")
# For pretty printing, replace sympy symbols with string representations
final_inequality_str = str(final_inequality)
final_inequality_str = final_inequality_str.replace('Le', '<=')
final_inequality_str = final_inequality_str.replace('alpha_i^D', 'alpha^D_i')
final_inequality_str = final_inequality_str.replace('(K*alpha^D)_i', '(K alpha^D)_i')
final_inequality_str = final_inequality_str.replace('(K*alpha^(D-i))_i', '(K alpha^(D-i))_i')
# Print the final formatted string
lhs_str = f"-{K_alpha_D_minus_i_i}".replace('(K*alpha^(D-i))_i', '(K alpha^(D-i))_i')
rhs_str = f"(1 + ({c1_val})*beta)*alpha^D_i - (1 + {c2_val}*beta)*(K alpha^D)_i + o(beta)"
print(f"{lhs_str} <= {rhs_str}")
