import math

# 1. Define constants and given parameters
# Loads
Gk = 1000  # kN (Permanent vertical load)
Qk_v = 1500  # kN (Variable vertical load)
Qk_h = 300  # kN (Variable horizontal load)

# Partial factors for actions (EC7 DA1-1)
gamma_G = 1.35
gamma_Q = 1.50

# Footing and column geometry
h_col = 2.0  # m (Height of column above ground)
D = 0.75  # m (Depth of footing base, thickness)
gamma_concrete = 24  # kN/m^3

# Soil properties
phi_k_deg = 35  # degrees (Characteristic angle of shearing resistance)
c_k = 0  # kPa (Characteristic effective cohesion)
gamma_soil = 20  # kN/m^3

# Partial factors for material properties (EC7 DA1-1)
gamma_phi = 1.0
gamma_c = 1.0
gamma_gamma = 1.0

# --- Iteratively determined footing size ---
# Through an iterative process, B=2.3m is found to be the optimal size
# where design action is approximately equal to design resistance.
B = 2.3  # m
L = B  # m (Square footing)

# 2. Calculate Design Actions (Loads)
Vd_superstructure = gamma_G * Gk + gamma_Q * Qk_v
Hd = gamma_Q * Qk_h

# Include self-weight of the footing
W_footing_k = B * L * D * gamma_concrete
Wd_footing = gamma_G * W_footing_k
Vd_total = Vd_superstructure + Wd_footing

# 3. Calculate Eccentricity and Effective Dimensions
M = Hd * (h_col + D)
e = M / Vd_total
B_prime = B - 2 * e
L_prime = L
A_prime = B_prime * L_prime

# 4. Calculate Design Bearing Resistance (R/A')
# Design soil parameters
phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi))
phi_d_rad = math.radians(phi_d_deg)
c_d = c_k / gamma_c
gamma_d = gamma_soil / gamma_gamma

# Effective overburden pressure at foundation base
q_prime = gamma_d * D

# Bearing capacity factors (EC7 Annex D)
Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
Nc = (Nq - 1) / math.tan(phi_d_rad) if math.tan(phi_d_rad) != 0 else 0
Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)

# Shape factors
sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
sgamma = 1 - 0.3 * (B_prime / L_prime)
sc = (sq * Nq - 1) / (Nq - 1) if (Nq - 1) != 0 else 1

# Inclination factors
m = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
# Vd_total is V in the formula, Hd is H
iq = (1 - Hd / Vd_total)**m
igamma = (1 - Hd / Vd_total)**(m + 1)
ic = iq - (1 - iq) / (Nc * math.tan(phi_d_rad)) if (Nc * math.tan(phi_d_rad)) != 0 else 0

# Design bearing resistance R/A'
# R/A' = c'd*Nc*sc*ic + q'*Nq*sq*iq + 0.5*γ'*B'*Nγ*sγ*iγ
term1 = c_d * Nc * sc * ic
term2 = q_prime * Nq * sq * iq
term3 = 0.5 * gamma_d * B_prime * Ngamma * sgamma * igamma
R_div_A = term1 + term2 + term3

# 5. Print the final result
print("The required ULS design bearing resistance is calculated using the formula:")
print("R/A' = c'd*Nc*sc*ic + q'*Nq*sq*iq + 0.5*γ'*B'*Nγ*sγ*iγ")
print("\nFor a footing of size B = {:.2f} m:".format(B))
print("\nSubstituting the calculated values into the resistance equation (since c'=0):")
print("R/A' = ({:.2f} * {:.2f} * {:.2f} * {:.2f}) + (0.5 * {:.2f} * {:.2f} * {:.2f} * {:.2f} * {:.2f})".format(
    q_prime, Nq, sq, iq, gamma_d, B_prime, Ngamma, sgamma, igamma))

# For comparison, let's also calculate the design bearing pressure
q_d = Vd_total / A_prime
print("\n--- Verification ---")
print(f"Design Bearing Pressure (Action), q_d = Vd / A' = {q_d:.2f} kN/m^2")
print(f"Design Bearing Resistance, R/A' = {R_div_A:.2f} kN/m^2")
print("\nSince the action and resistance are almost equal, the footing size is optimal.")
print("\nThe required ULS design bearing resistance is:")
print(f"{R_div_A:.1f}")