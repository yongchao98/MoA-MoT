import math

# 1. Define Parameters and Assumptions
# Given characteristic values
G_k = 500.0  # kN (Permanent vertical load)
Q_v_k = 1500.0 # kN (Variable vertical load)
Q_h_k = 120.0  # kN (Variable horizontal load)
zeta = 0.7     # Combination factor for horizontal load

# Foundation and soil properties
B = 2.39  # m (Foundation width)
L = 2.39  # m (Foundation length)
D = 1.0   # m (Foundation depth)
gamma_soil = 20.0  # kN/m^3 (Soil unit weight)
c_k_prime = 10.0   # kPa or kN/m^2 (Characteristic cohesion)
phi_k_prime_deg = 28.0 # degrees (Characteristic friction angle)

# Assumed partial factors (EC7, DA1-Combination 1 for actions, DA1-Combination 2 for materials)
gamma_G = 1.35
gamma_Q = 1.5
gamma_c_prime = 1.25
gamma_phi_prime = 1.25

print("--- Step 1: Design Material Properties ---")
# 2. Calculate Design Material Properties
c_d_prime = c_k_prime / gamma_c_prime
phi_k_prime_rad = math.radians(phi_k_prime_deg)
phi_d_prime_rad = math.atan(math.tan(phi_k_prime_rad) / gamma_phi_prime)
phi_d_prime_deg = math.degrees(phi_d_prime_rad)
print(f"Design cohesion (c'_d): {c_d_prime:.2f} kPa")
print(f"Design friction angle (φ'_d): {phi_d_prime_deg:.2f} degrees")

print("\n--- Step 2: Design Loads for Inclination Factors ---")
# 3. Calculate Design Loads
V_d = gamma_G * G_k + gamma_Q * Q_v_k
# The horizontal load is a variable action, so it gets gamma_Q.
# It's an accompanying action, so it's multiplied by the combination factor zeta.
H_d = gamma_Q * zeta * Q_h_k
print(f"Design vertical load (V_d): {V_d:.2f} kN")
print(f"Design horizontal load (H_d): {H_d:.2f} kN")

print("\n--- Step 3: Effective Foundation Area ---")
B_prime = B
L_prime = L
A_prime = B_prime * L_prime
print(f"Effective area (A'): {A_prime:.4f} m^2")

print("\n--- Step 4: Bearing Capacity Factors ---")
# 4. Calculate Bearing Capacity Factors (based on φ'_d)
N_q = math.exp(math.pi * math.tan(phi_d_prime_rad)) * math.pow(math.tan(math.pi/4 + phi_d_prime_rad/2), 2)
N_c = (N_q - 1) / math.tan(phi_d_prime_rad) if math.tan(phi_d_prime_rad) > 0 else (2 + math.pi)
N_gamma = 2 * (N_q - 1) * math.tan(phi_d_prime_rad)
print(f"N_q = {N_q:.3f}")
print(f"N_c = {N_c:.3f}")
print(f"N_γ = {N_gamma:.3f}")

print("\n--- Step 5: Shape, Depth, and Inclination Factors ---")
# Shape factors
s_q = 1 + (B_prime / L_prime) * math.sin(phi_d_prime_rad)
s_gamma = 1 - 0.3 * (B_prime / L_prime)
s_c = (s_q * N_q - 1) / (N_q - 1)
print(f"Shape factors: s_q = {s_q:.3f}, s_c = {s_c:.3f}, s_γ = {s_gamma:.3f}")

# Depth factors
d_q = 1 + 2 * math.tan(phi_d_prime_rad) * math.pow(1 - math.sin(phi_d_prime_rad), 2) * (D / B)
d_gamma = 1.0
d_c = d_q - (1 - d_q) / (N_c * math.tan(phi_d_prime_rad))
print(f"Depth factors: d_q = {d_q:.3f}, d_c = {d_c:.3f}, d_γ = {d_gamma:.3f}")

# Inclination factors
m = (2 + B_prime / L_prime) / (1 + B_prime / L_prime)
# As per problem, omit A'c'cot(phi) term in denominator
i_q = math.pow(1 - H_d / V_d, m)
i_gamma = math.pow(1 - H_d / V_d, m + 1)
i_c = i_q - (1 - i_q) / (N_c * math.tan(phi_d_prime_rad))
print(f"Inclination factors: i_q = {i_q:.3f}, i_c = {i_c:.3f}, i_γ = {i_gamma:.3f}")

print("\n--- Step 6: Calculate Design Bearing Resistance (R_d) ---")
# Overburden pressure at foundation base
q = gamma_soil * D
print(f"Overburden pressure (q): {q:.2f} kPa")

# Calculate the three terms of the bearing capacity equation
term_c = c_d_prime * N_c * s_c * d_c * i_c
term_q = q * N_q * s_q * d_q * i_q
term_gamma = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * d_gamma * i_gamma

# Calculate total design resistance
R_d = A_prime * (term_c + term_q + term_gamma)

print("\nThe design resistance R_d is calculated as:")
print(f"R_d = A' * (c'_d*N_c*s_c*d_c*i_c + q*N_q*s_q*d_q*i_q + 0.5*γ*B'*N_γ*s_γ*d_γ*i_γ)")
print(f"R_d = {A_prime:.4f} * ({c_d_prime:.2f}*{N_c:.3f}*{s_c:.3f}*{d_c:.3f}*{i_c:.3f} + "
      f"{q:.2f}*{N_q:.3f}*{s_q:.3f}*{d_q:.3f}*{i_q:.3f} + "
      f"0.5*{gamma_soil:.2f}*{B_prime:.2f}*{N_gamma:.3f}*{s_gamma:.3f}*{d_gamma:.3f}*{i_gamma:.3f})")
print(f"R_d = {A_prime:.4f} * ({term_c:.2f} + {term_q:.2f} + {term_gamma:.2f})")
print(f"R_d = {A_prime:.4f} * ({term_c + term_q + term_gamma:.2f})")
print(f"\nFinal Design Resistance (R_d) = {R_d:.1f} kN")