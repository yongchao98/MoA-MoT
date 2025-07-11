import sympy

# Define the symbolic variables
V, epsilon, pi, L, a, b = sympy.symbols('V varepsilon pi L a b')
sigma_0 = sympy.symbols('sigma_0')

# Formulas from Answer Choice B
q_v_B = (-4 * V * epsilon * pi * L) / (1 - a**2/b**2)
q_s_a_B = (2 * pi * L * V * epsilon) / (1 - a**2/b**2)
q_s_b_B = (-4 * pi * L * V * epsilon * a**2) / (b**2 * (1 - a**2/b**2))

# Print the expressions
print("Formulas from the selected answer choice:")
print(f"Total volume charge = q_v = {sympy.pretty(q_v_B)}")
print(f"Total surface charge on inner electrode = q_s(r = a) = {sympy.pretty(q_s_a_B)}")
print(f"Total surface charge on outer electrode = q_s(r = b) = {sympy.pretty(q_s_b_B)}")