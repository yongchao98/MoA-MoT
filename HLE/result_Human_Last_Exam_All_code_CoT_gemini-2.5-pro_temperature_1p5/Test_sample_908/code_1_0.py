import sympy as sp

# Define the symbols used in the equations
V, epsilon, L, a, b = sp.symbols('V epsilon L a b')
pi = sp.pi

# --- Formulas from Answer Choice B ---

# Common denominator term
denominator = 1 - (a**2 / b**2)

# Total volume charge from option B
q_v_B = -4 * V * epsilon * pi * L / denominator

# Total surface charge on the inner electrode (r=a) from option B
q_s_a_B = 2 * pi * L * V * epsilon / denominator

# Total surface charge on the outer electrode (r=b) from option B
q_s_b_B = -4 * pi * L * V * epsilon * a**2 / (b**2 * denominator)

# Print the final equations
print("The equations corresponding to Answer Choice B are:")
print("Total volume charge (q_v):")
# The following line prints the equation for q_v
print(f"q_v = {sp.printing.latex(q_v_B)}")


print("\nTotal surface charge on inner electrode (q_s(a)):")
# The following line prints the equation for q_s(a)
print(f"q_s(a) = {sp.printing.latex(q_s_a_B)}")


print("\nTotal surface charge on outer electrode (q_s(b)):")
# The following line prints the equation for q_s(b)
print(f"q_s(b) = {sp.printing.latex(q_s_b_B)}")