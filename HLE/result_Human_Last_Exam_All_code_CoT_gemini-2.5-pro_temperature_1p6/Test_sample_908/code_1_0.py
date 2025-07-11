import sympy as sp

# Define the symbols
V, epsilon, pi, L, a, b = sp.symbols('V epsilon pi L a b', positive=True)

# Formulas from Option B
q_v_expr = (-4 * V * epsilon * pi * L) / (1 - a**2 / b**2)
qs_a_expr = (2 * pi * L * V * epsilon) / (1 - a**2 / b**2)
qs_b_expr = (-4 * pi * L * V * epsilon * a**2) / (b**2 * (1 - a**2/b**2))

# Print the results
print("The solution based on the provided answer choices is:")
print("Total volume charge = q_v = ", end="")
sp.pretty_print(q_v_expr)
print("\nTotal surface charge on inner electrode = q_s(r = a) = ", end="")
sp.pretty_print(qs_a_expr)
print("\nTotal surface charge on outer electrode = q_s(r = b) = ", end="")
sp.pretty_print(qs_b_expr)
