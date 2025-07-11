import math

# Step 1: Define Input Parameters
# Loads
G = 500.0  # Permanent vertical load (kN)
Q_v = 1500.0 # Variable vertical load (kN)
Q_h = 120.0  # Variable horizontal load (kN)

# Factors
gamma_G = 1.35 # Partial factor for permanent actions
gamma_Q = 1.50 # Partial factor for variable actions
zeta = 0.7   # Combination factor for independent variable loads

# Foundation and Soil Properties
B = 2.39  # Foundation width (m)
L = 2.39  # Foundation length (m)
D = 1.0   # Foundation depth (m)
c_k = 10.0  # Characteristic cohesion (kPa or kN/m^2)
phi_k_deg = 28.0 # Characteristic friction angle (degrees)
gamma_soil = 20.0 # Soil unit weight (kN/m^3)

# Step 2: Calculate Design Actions (Loads)
V_d = gamma_G * G + gamma_Q * Q_v
H_d = gamma_Q * zeta * Q_h

# Step 3: Determine Effective Foundation Area
# Moment M_d is caused by H_d acting at the top of the foundation (lever arm = D)
M_d = H_d * D
# Eccentricity e
e = M_d / V_d
# Effective dimensions
B_prime = B - 2 * e
L_prime = L # Eccentricity is along the width B
A_prime = B_prime * L_prime

# Step 4: Determine Design Soil Parameters (DA1-1)
gamma_c_prime = 1.0
gamma_phi_prime = 1.0
c_d = c_k / gamma_c_prime
phi_d_rad = math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi_prime)
phi_d_deg = math.degrees(phi_d_rad)

# Step 5: Calculate Bearing Capacity Factors (for phi_d)
tan_phi_d = math.tan(phi_d_rad)
tan_sq_45_plus_phi_half = math.tan(math.radians(45) + phi_d_rad / 2)**2

N_q = math.exp(math.pi * tan_phi_d) * tan_sq_45_plus_phi_half
N_c = (N_q - 1) / tan_phi_d if tan_phi_d != 0 else 0
N_gamma = 2 * (N_q - 1) * tan_phi_d # Vesic's formula

# Step 6: Calculate Shape and Inclination Factors
# Overburden pressure at foundation base
q_prime = gamma_soil * D

# Shape factors
s_q = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
s_gamma = 1 - 0.3 * (B_prime / L_prime)
s_c = (s_q * N_q - 1) / (N_q - 1) if (N_q - 1) != 0 else 1

# Inclination factors
# Exponent m
m = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
# As per instruction, omit A'c'cot(phi) term from denominator
i_q = (1 - H_d / V_d)**m
i_gamma = (1 - H_d / V_d)**(m + 1)
i_c = i_q - (1 - i_q) / (N_c * tan_phi_d) if (N_c * tan_phi_d) != 0 else i_q

# Step 7: Calculate Design Bearing Resistance
# Calculate each term of the bearing capacity equation
term_c = c_d * N_c * s_c * i_c
term_q = q_prime * N_q * s_q * i_q
term_gamma = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * i_gamma

# Ultimate bearing pressure (q_ult)
q_ult = term_c + term_q + term_gamma

# Total design resistance (R_d)
R_d = q_ult * A_prime

# Step 8: Output the Result
print("Calculation of Design Resistance (R_d)\n")
print("1. Design Loads:")
print(f"V_d = {gamma_G} * {G} + {gamma_Q} * {Q_v} = {V_d:.2f} kN")
print(f"H_d = {gamma_Q} * {zeta} * {Q_h} = {H_d:.2f} kN\n")

print("2. Effective Dimensions:")
print(f"e = M_d / V_d = ({H_d:.2f} * {D}) / {V_d:.2f} = {e:.4f} m")
print(f"B' = B - 2*e = {B} - 2*{e:.4f} = {B_prime:.4f} m")
print(f"A' = B' * L' = {B_prime:.4f} * {L_prime} = {A_prime:.4f} m^2\n")

print("3. Factors:")
print(f"Bearing Capacity Factors (φ'_d = {phi_d_deg:.2f}°):")
print(f"N_c = {N_c:.2f}, N_q = {N_q:.2f}, N_γ = {N_gamma:.2f}")
print("Shape Factors:")
print(f"s_c = {s_c:.3f}, s_q = {s_q:.3f}, s_γ = {s_gamma:.3f}")
print("Inclination Factors:")
print(f"i_c = {i_c:.3f}, i_q = {i_q:.3f}, i_γ = {i_gamma:.3f}\n")

print("4. Design Resistance Calculation:")
print("R_d = A' * (c'_d*N_c*s_c*i_c + q'*N_q*s_q*i_q + 0.5*γ'*B'*N_γ*s_γ*i_γ)")
print(f"R_d = {A_prime:.3f} * ({c_d:.1f}*{N_c:.2f}*{s_c:.3f}*{i_c:.3f} + "
      f"{q_prime:.1f}*{N_q:.2f}*{s_q:.3f}*{i_q:.3f} + "
      f"0.5*{gamma_soil:.1f}*{B_prime:.3f}*{N_gamma:.2f}*{s_gamma:.3f}*{i_gamma:.3f})")
print(f"R_d = {A_prime:.3f} * ({term_c:.1f} + {term_q:.1f} + {term_gamma:.1f})")
print(f"R_d = {A_prime:.3f} * {q_ult:.1f}")
print(f"R_d = {R_d:.1f} kN")

# Final answer in the required format
final_answer = R_d
print(f"\n<<< {final_answer:.1f} >>>")