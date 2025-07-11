import sympy

# Define the symbolic variables
V, epsilon, pi, L, a, b = sympy.symbols('V varepsilon pi L a b')

# Expressions from Option B
q_v_expr_B = (-4 * V * epsilon * pi * L) / (1 - a**2/b**2)
q_s_a_expr_B = (2 * pi * L * V * epsilon) / (1 - a**2/b**2)
q_s_b_expr_B = (-4 * pi * L * V * epsilon * a**2) / (b**2 * (1 - a**2/b**2))

# Use pretty print for better visualization
print("The solution based on the provided options is:")
print("-" * 20)
print("Total volume charge = q_v = " + sympy.pretty(q_v_expr_B, use_unicode=False))
print("")
print("Total surface charge on inner electrode = q_s(r = a) = " + sympy.pretty(q_s_a_expr_B, use_unicode=False))
print("")
print("Total surface charge on outer electrode = q_s(r = b) = " + sympy.pretty(q_s_b_expr_B, use_unicode=False))
