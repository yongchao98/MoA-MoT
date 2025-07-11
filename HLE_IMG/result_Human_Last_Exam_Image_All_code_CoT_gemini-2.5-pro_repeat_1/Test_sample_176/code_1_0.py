import math

# Step 1: Define Inputs
G_k = 500.0  # Permanent vertical load [kN]
Q_v_k = 1500.0 # Variable vertical load [kN]
Q_h_k = 120.0  # Variable horizontal load [kN]
zeta = 0.7     # Combination factor
B = 2.39       # Foundation width [m]
L = 2.39       # Foundation length [m]
D = 1.0        # Foundation depth [m]
h = 1.0        # Assumed foundation thickness [m]
gamma_c = 24.0 # Concrete unit weight [kN/m^3]
c_k = 10.0     # Characteristic cohesion [kPa or kN/m^2]
phi_k = 28.0   # Characteristic friction angle [degrees]
gamma_s = 20.0 # Soil unit weight [kN/m^3]

# Partial factors for DA1-1
gamma_G = 1.35
gamma_Q = 1.5
gamma_c_prime = 1.0
gamma_phi_prime = 1.0
gamma_R = 1.0

print("Step 1: Given Parameters")
print(f"Permanent load G_k = {G_k} kN")
print(f"Variable vertical load Q_v_k = {Q_v_k} kN")
print(f"Variable horizontal load Q_h_k = {Q_h_k} kN")
print(f"Foundation width B = {B} m, Length L = {L} m, Depth D = {h} m")
print(f"Soil properties: c'_k = {c_k} kPa, phi'_k = {phi_k} degrees, gamma_s = {gamma_s} kN/m^3\n")

# Step 2: Calculate Design Actions (Loads)
W_f = B * L * h * gamma_c
G_total_k = G_k + W_f
V_d = gamma_G * G_total_k + gamma_Q * Q_v_k
H_d = gamma_Q * zeta * Q_h_k

print("Step 2: Design Actions (Loads)")
print(f"Foundation self-weight W_f = {B} * {L} * {h} * {gamma_c} = {W_f:.2f} kN")
print(f"Total design vertical load V_d = {gamma_G} * ({G_k} + {W_f:.2f}) + {gamma_Q} * {Q_v_k} = {V_d:.2f} kN")
print(f"Design horizontal load H_d = {gamma_Q} * {zeta} * {Q_h_k} = {H_d:.2f} kN\n")

# Step 3: Calculate Effective Foundation Area
M_d = H_d * h
e_B = M_d / V_d
B_prime = B - 2 * e_B
L_prime = L
A_prime = B_prime * L_prime

print("Step 3: Effective Foundation Dimensions")
print(f"Moment at base M_d = {H_d:.2f} kN * {h} m = {M_d:.2f} kNm")
print(f"Eccentricity e_B = {M_d:.2f} / {V_d:.2f} = {e_B:.4f} m")
print(f"Effective width B' = {B} - 2 * {e_B:.4f} = {B_prime:.4f} m")
print(f"Effective length L' = {L_prime:.4f} m")
print(f"Effective area A' = {B_prime:.4f} * {L_prime:.4f} = {A_prime:.4f} m^2\n")

# Step 4: Determine Design Soil Parameters
c_d = c_k / gamma_c_prime
phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k)) / gamma_phi_prime))
phi_d_rad = math.radians(phi_d_deg)

print("Step 4: Design Soil Parameters (DA1-1, gamma_M=1.0)")
print(f"Design cohesion c'_d = {c_k} / {gamma_c_prime} = {c_d:.2f} kPa")
print(f"Design friction angle phi'_d = {phi_d_deg:.2f} degrees\n")

# Step 5: Calculate Bearing Capacity Factors
N_q = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
N_c = (N_q - 1) / math.tan(phi_d_rad)
N_gamma = 2 * (N_q - 1) * math.tan(phi_d_rad)

print("Step 5: Bearing Capacity Factors")
print(f"N_q = {N_q:.3f}")
print(f"N_c = {N_c:.3f}")
print(f"N_gamma = {N_gamma:.3f}\n")

# Step 6: Calculate Shape Factors
s_q = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
s_gamma = 1 - 0.3 * (B_prime / L_prime)
s_c = (s_q * N_q - 1) / (N_q - 1)

print("Step 6: Shape Factors")
print(f"s_q = 1 + ({B_prime:.3f}/{L_prime:.3f}) * sin({phi_d_deg:.1f}) = {s_q:.3f}")
print(f"s_gamma = 1 - 0.3 * ({B_prime:.3f}/{L_prime:.3f}) = {s_gamma:.3f}")
print(f"s_c = ({s_q:.3f} * {N_q:.3f} - 1) / ({N_q:.3f} - 1) = {s_c:.3f}\n")

# Step 7: Calculate Load Inclination Factors
m = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
i_q = (1 - H_d / V_d)**m
i_c = i_q
i_gamma = (1 - H_d / V_d)**(m + 1)

print("Step 7: Load Inclination Factors")
print(f"Exponent m = (2 + {B_prime:.3f}/{L_prime:.3f}) / (1 + {B_prime:.3f}/{L_prime:.3f}) = {m:.3f}")
print(f"i_q = (1 - {H_d:.2f}/{V_d:.2f})^{m:.3f} = {i_q:.3f}")
print(f"i_c = i_q = {i_c:.3f}")
print(f"i_gamma = (1 - {H_d:.2f}/{V_d:.2f})^({m:.3f}+1) = {i_gamma:.3f}\n")

# Step 8: Calculate Design Bearing Resistance
q_prime = gamma_s * D
term_c = c_d * N_c * s_c * i_c
term_q = q_prime * N_q * s_q * i_q
term_gamma = 0.5 * gamma_s * B_prime * N_gamma * s_gamma * i_gamma
R_d = A_prime * (term_c + term_q + term_gamma)

print("Step 8: Design Bearing Resistance Calculation")
print(f"Surcharge q' = {gamma_s} * {D} = {q_prime:.2f} kPa")
print("\nFinal Bearing Resistance Equation:")
print(f"R_d = A' * (c'*N_c*s_c*i_c + q'*N_q*s_q*i_q + 0.5*gamma*B'*N_gamma*s_gamma*i_gamma)")
print(f"R_d = {A_prime:.3f} * (({c_d:.1f}*{N_c:.2f}*{s_c:.3f}*{i_c:.3f}) + ({q_prime:.1f}*{N_q:.2f}*{s_q:.3f}*{i_q:.3f}) + (0.5*{gamma_s:.1f}*{B_prime:.3f}*{N_gamma:.2f}*{s_gamma:.3f}*{i_gamma:.3f}))")
print(f"R_d = {A_prime:.3f} * ({term_c:.2f} + {term_q:.2f} + {term_gamma:.2f})")
print(f"R_d = {A_prime:.3f} * ({term_c + term_q + term_gamma:.2f})")
print(f"Total Design Resistance R_d = {R_d:.1f} kN")

final_answer = round(R_d, 1)
# print(f"\nFinal Answer: {final_answer}")