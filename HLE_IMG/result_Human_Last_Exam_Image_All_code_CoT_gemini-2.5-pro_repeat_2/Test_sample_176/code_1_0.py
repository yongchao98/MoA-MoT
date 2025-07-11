import math

# Step 1: Define Inputs
G_k = 500.0  # Permanent vertical load [kN]
Qv_k = 1500.0 # Variable vertical load [kN]
Qh_k = 120.0  # Variable horizontal load [kN]
zeta = 0.7    # Combination factor

B = 2.39      # Foundation width [m]
L = 2.39      # Foundation length [m]
D = 1.0       # Foundation depth [m]

gamma_c = 24.0 # Unit weight of concrete [kN/m^3]
gamma_s = 20.0 # Unit weight of soil [kN/m^3]
c_k = 10.0    # Characteristic cohesion [kPa or kN/m^2]
phi_k_deg = 28.0 # Characteristic friction angle [degrees]

# Step 2: Set Partial Factors (Eurocode 7, DA1, Combination 1)
gamma_G = 1.35
gamma_Q = 1.50
gamma_c_prime = 1.0
gamma_phi_prime = 1.0
gamma_gamma = 1.0
gamma_R_v = 1.0 # Partial factor for bearing resistance

# Step 3: Calculate Design Actions
# Assume foundation thickness h = D = 1.0 m based on diagram
h = D
W_fk = B * L * h * gamma_c
W_fd = gamma_G * W_fk

Vd_applied = gamma_G * G_k + gamma_Q * Qv_k
Vd_total = Vd_applied + W_fd

Hd = gamma_Q * zeta * Qh_k
# Lever arm is the depth D, as load is applied at the top of the foundation (ground level)
Md = Hd * D

# Step 4: Determine Effective Dimensions
e_B = Md / Vd_total
B_prime = B - 2 * e_B
L_prime = L # No eccentricity in L direction
A_prime = B_prime * L_prime

# Step 5: Calculate Design Soil Parameters
c_d = c_k / gamma_c_prime
phi_k_rad = math.radians(phi_k_deg)
phi_d_rad = math.atan(math.tan(phi_k_rad) / gamma_phi_prime)
phi_d_deg = math.degrees(phi_d_rad)
gamma_d = gamma_s / gamma_gamma
q_prime = gamma_d * D

# Step 6: Compute Bearing Capacity Factors (for phi_d)
Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
Nc = (Nq - 1) / math.tan(phi_d_rad) if math.tan(phi_d_rad) > 0 else (2 + math.pi)
N_gamma = 2 * (Nq - 1) * math.tan(phi_d_rad)

# Step 7: Compute Shape Factors
ratio_B_L = B_prime / L_prime
s_q = 1 + ratio_B_L * math.sin(phi_d_rad)
s_gamma = 1 - 0.3 * ratio_B_L
s_c = (s_q * Nq - 1) / (Nq - 1)

# Step 8: Compute Inclination Factors
# Exponent m for load inclined parallel to B'
m = (2 + ratio_B_L) / (1 + ratio_B_L)
# Per prompt, omit A'c'cot(phi) term from denominator
iq = (1 - Hd / Vd_total)**m
i_gamma = (1 - Hd / Vd_total)**(m + 1)
ic = (iq * Nq - 1) / (Nq - 1)

# Step 9: Calculate Design Resistance (R_d)
term_c = c_d * Nc * s_c * ic
term_q = q_prime * Nq * s_q * iq
term_gamma = 0.5 * gamma_d * B_prime * N_gamma * s_gamma * i_gamma

q_ult = term_c + term_q + term_gamma
R_d = q_ult * A_prime / gamma_R_v

# Print the final equation and result
print("Calculation of Design Bearing Resistance (R_d)\n")
print("R_d = (c'_d * N_c * s_c * i_c + q' * N_q * s_q * i_q + 0.5 * γ'_d * B' * N_γ * s_γ * i_γ) * A' / γ_R,v\n")

print(f"R_d = ({c_d:.2f} * {Nc:.2f} * {s_c:.3f} * {ic:.3f} + "
      f"{q_prime:.2f} * {Nq:.2f} * {s_q:.3f} * {iq:.3f} + "
      f"0.5 * {gamma_d:.2f} * {B_prime:.3f} * {N_gamma:.2f} * {s_gamma:.3f} * {i_gamma:.3f}) * "
      f"{A_prime:.3f} / {gamma_R_v:.1f}\n")

print(f"R_d = ({term_c:.2f} + {term_q:.2f} + {term_gamma:.2f}) * {A_prime:.3f}\n")
print(f"R_d = {q_ult:.2f} * {A_prime:.3f}\n")
print(f"The design resistance R_d is: {R_d:.1f} kN")
