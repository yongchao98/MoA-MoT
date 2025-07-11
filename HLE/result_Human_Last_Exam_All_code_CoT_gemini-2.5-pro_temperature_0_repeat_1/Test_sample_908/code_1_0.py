import sympy as sp

# Define the symbols
V, epsilon, pi, L, a, b = sp.symbols('V epsilon pi L a b', real=True, positive=True)

# Formulas from option B
q_v_B = -4 * V * epsilon * pi * L / (1 - a**2/b**2)
q_s_a_B = 2 * pi * L * V * epsilon / (1 - a**2/b**2)
q_s_b_B = -4 * pi * L * V * epsilon * a**2 / (b**2 * (1 - a**2/b**2))

# Print the expressions
print("Total volume charge = q_v = {}".format(sp.pretty(q_v_B, use_unicode=False)))
print("Total surface charge on inner electrode = q_s(r = a) = {}".format(sp.pretty(q_s_a_B, use_unicode=False)))
print("Total surface charge on outer electrode = q_s(r = b) = {}".format(sp.pretty(q_s_b_B, use_unicode=False)))
