import sympy

# Define variables for the equation
alpha_i, K_alpha_i, K_alpha_loo_i, beta = sympy.symbols("alpha_i^D (K*alpha^D)_i (K*alpha^{D-i})_i beta")
c1, c2 = -1, 0

# Construct the RHS of the inequality
rhs = (1 + c1 * beta) * alpha_i - (1 + c2 * beta) * K_alpha_i

# Display the final inequality
print("The determined coefficients are c1 = -1 and c2 = 0.")
print("The extended Jaakola-Haussler bound is:")

# Using sympy.pretty_print for a clearer mathematical output
inequality = sympy.Le(-K_alpha_loo_i, rhs)
sympy.pretty_print(inequality, use_unicode=True)

# Also printing each number as requested
print("\nFinal equation with explicit numbers:")
print(f"-{K_alpha_loo_i} <= (1 + ({c1})*β)*{alpha_i} - (1 + ({c2})*β)*{K_alpha_i}")
