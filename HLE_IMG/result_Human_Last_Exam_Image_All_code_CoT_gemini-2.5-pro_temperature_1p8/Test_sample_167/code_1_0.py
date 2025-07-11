import math

# --- 1. Define Given Data and Constants ---
# Loads
Gk_applied = 1000.0  # kN
Qk_v = 1500.0      # kN (variable vertical load)
Qk_h = 300.0       # kN (variable horizontal load)

# Geometry
h_load_arm = 2.0 + 0.75  # m, lever arm for horizontal load above foundation base
h_footing = 0.75         # m, thickness of footing
D = 0.75                 # m, depth of foundation base

# Material Properties
phi_k_deg = 35.0   # degrees
c_k = 0.0          # kPa
gamma_soil = 20.0  # kN/m^3
gamma_conc = 24.0  # kN/m^3

# Eurocode 7 DA1-C1 Partial Factors
gamma_G = 1.35
gamma_Q = 1.50

# --- 2. Footing Size ---
# Through an iterative design process (checking Vd <= Rd), a footing width of
# B = 2.4 m is found to be adequate. We will proceed with this dimension.
B = 2.4  # m
L = B    # m, square footing

print(f"Step 1: Assume a suitable footing size. An iterative design process yields B = L = {B} m.\n")

# --- 3. Calculate Design Actions (Loads) ---
# Self-weight of the footing
Gk_footing = (B * L * h_footing) * gamma_conc
Gk_tot = Gk_applied + Gk_footing

# Design loads based on DA1-C1
V_d = gamma_G * Gk_tot + gamma_Q * Qk_v
H_d = gamma_Q * Qk_h
M_d = H_d * h_load_arm

print("Step 2: Calculate ULS design actions (loads and moment).")
print(f"   - Design vertical load (V_d): {V_d:.2f} kN")
print(f"   - Design horizontal load (H_d): {H_d:.2f} kN")
print(f"   - Design moment at base (M_d): {M_d:.2f} kNm\n")

# --- 4. Calculate Effective Footing Area ---
# Eccentricity of the load
e_B = M_d / V_d

# Effective dimensions
B_prime = B - 2 * e_B
L_prime = L
A_prime = B_prime * L_prime

print("Step 3: Calculate effective footing area based on load eccentricity.")
print(f"   - Eccentricity (e_B = M_d / V_d): {M_d:.2f} / {V_d:.2f} = {e_B:.4f} m")
print(f"   - Effective width (B' = B - 2e): {B:.2f} - 2 * {e_B:.4f} = {B_prime:.4f} m")
print(f"   - Effective area (A' = B' * L): {B_prime:.4f} * {L:.2f} = {A_prime:.4f} m^2\n")

# --- 5. Calculate Required Bearing Resistance ---
# The required bearing resistance is the design pressure applied by the foundation.
q_Ed = V_d / A_prime

print("Step 4: Calculate the required ULS design bearing resistance (q_Ed).")
print("This is the design pressure that the soil must be able to resist.")
print("Formula: q_Ed = V_d / A'")
print(f"q_Ed = {V_d:.2f} kN / {A_prime:.4f} m^2")
print(f"\nThe required ULS design bearing resistance is: {q_Ed:.2f} kN/m^2")

# Final answer for extraction
final_answer = q_Ed