import sympy as sp

# Define the symbols used in the problem
V, epsilon, pi, L, a, b = sp.symbols('V varepsilon pi L a b', real=True, positive=True)

# From the derivation, the constant C is:
C = (2 * V * a**2 * b**2) / (b**2 - a**2)

# 1. Calculate the total volume charge (q_v)
# The volume charge density is rho_v = -2 * epsilon * C / r**4
# The total volume charge is the integral of rho_v over the volume
# q_v = Integral(-2*epsilon*C/r**4 * 2*pi*r*L, (r, a, b))
# This simplifies to:
q_v = -4 * pi * L * epsilon * V

# 2. Calculate the total surface charge on the inner electrode (q_s_a)
# The surface charge density is rho_s_a = epsilon * E(a) = epsilon * C / a**3
# The total surface charge is rho_s_a * Area = (epsilon * C / a**3) * (2 * pi * a * L)
q_s_a = (2 * pi * L * epsilon * C) / a**2
# Substitute C and simplify
q_s_a_simplified = sp.simplify(q_s_a)

# 3. Calculate the total surface charge on the outer electrode (q_s_b)
# The surface charge density is rho_s_b = -epsilon * E(b) = -epsilon * C / b**3
# The total surface charge is rho_s_b * Area = (-epsilon * C / b**3) * (2 * pi * b * L)
q_s_b = (-2 * pi * L * epsilon * C) / b**2
# Substitute C and simplify
q_s_b_simplified = sp.simplify(q_s_b)

# Print the derived expressions
print("Derived correct expressions:")
print("-" * 30)

print("Total volume charge (q_v):")
# The final equation has the number -4
print(f"q_v = {sp.pretty(q_v, use_unicode=False)}")
print("\n")

print("Total surface charge on inner electrode q_s(r=a):")
# The final equation has the number 4
# We format it to match the options
q_s_a_formatted = (4 * pi * L * V * epsilon) / (1 - a**2/b**2)
print(f"q_s(r=a) = {sp.pretty(q_s_a_formatted, use_unicode=False)}")
print("\n")

print("Total surface charge on outer electrode q_s(r=b):")
# The final equation has the number -4
# We format it to match the options
q_s_b_formatted = (-4 * pi * L * V * epsilon * a**2) / (b**2 * (1 - a**2/b**2))
print(f"q_s(r=b) = {sp.pretty(q_s_b_formatted, use_unicode=False)}")
print("-" * 30)
print("\nComparing these results with the answer choices, Option D has the correct expressions for the surface charges.")
