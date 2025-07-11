import math

# Step 1: Define Inputs
# Given parameters from the problem description
G_k = 500  # Permanent vertical load (kN)
Qv_k = 1500  # Variable vertical load (kN)
Qh_k = 120  # Variable horizontal load (kN)
zeta = 0.7  # Combination factor
B = 2.39  # Foundation width (m)
L = 2.39  # Foundation length (m)
Df = 1.0  # Foundation depth (m)
c_k = 10  # Characteristic effective cohesion (kPa or kN/m^2)
phi_k_deg = 28  # Characteristic effective friction angle (degrees)
gamma_soil = 20  # Soil unit weight (kN/m^3)

# Assumed parameters based on common practice and problem context
h = 1.0 # Assume foundation thickness is equal to its depth (m)
gamma_G = 1.35  # Partial factor for permanent action
gamma_Q = 1.50  # Partial factor for variable action
gamma_c = 1.0  # Material factor for cohesion
gamma_phi = 1.0 # Material factor for friction angle
gamma_gamma = 1.0 # Material factor for soil weight

# Step 2: Calculate Design Actions
V_d = gamma_G * G_k + gamma_Q * Qv_k
H_d = gamma_Q * zeta * Qh_k

print("--- Design Actions (Load Combination 1) ---")
print(f"Design vertical load V_d = {gamma_G} * {G_k} + {gamma_Q} * {Qv_k} = {V_d:.2f} kN")
print(f"Design horizontal load H_d = {gamma_Q} * {zeta} * {Qh_k} = {H_d:.2f} kN")
print("-" * 20)

# Step 3: Determine Effective Foundation Dimensions
M_d = H_d * h
e_B = M_d / V_d
B_prime = B - 2 * e_B
L_prime = L  # No eccentricity in L direction
A_prime = B_prime * L_prime

print("--- Effective Foundation Dimensions ---")
print(f"Moment at base M_d = {H_d:.2f} kN * {h:.2f} m = {M_d:.2f} kNm")
print(f"Eccentricity e_B = {M_d:.2f} kNm / {V_d:.2f} kN = {e_B:.4f} m")
print(f"Effective width B' = {B:.2f} m - 2 * {e_B:.4f} m = {B_prime:.4f} m")
print(f"Effective length L' = {L_prime:.2f} m")
print(f"Effective area A' = {B_prime:.4f} m * {L_prime:.2f} m = {A_prime:.4f} m^2")
print("-" * 20)

# Step 4: Calculate Design Soil Parameters
c_d = c_k / gamma_c
phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi))
phi_d_rad = math.radians(phi_d_deg)
gamma_d = gamma_soil / gamma_gamma

# Step 5: Calculate Bearing Capacity Factors
N_q = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.pi/4 + phi_d_rad/2))**2
N_c = (N_q - 1) / math.tan(phi_d_rad)
N_gamma = 2 * (N_q - 1) * math.tan(phi_d_rad)

print("--- Bearing Capacity Factors (for phi'_d = {:.1f} deg) ---".format(phi_d_deg))
print(f"N_q = {N_q:.3f}")
print(f"N_c = {N_c:.3f}")
print(f"N_gamma = {N_gamma:.3f}")
print("-" * 20)

# Surcharge pressure q'
q_prime = gamma_d * Df
print(f"Surcharge at foundation base q' = {gamma_d:.2f} kN/m^3 * {Df:.2f} m = {q_prime:.2f} kPa")
print("-" * 20)


# Step 6: Calculate Shape Factors
B_prime_over_L_prime = B_prime / L_prime
s_q = 1 + (B_prime_over_L_prime) * math.sin(phi_d_rad)
s_gamma = 1 - 0.3 * (B_prime_over_L_prime)
s_c = (s_q * N_q - 1) / (N_q - 1)

print("--- Shape Factors ---")
print(f"s_q = 1 + ({B_prime:.3f}/{L_prime:.3f}) * sin({phi_d_deg:.1f}) = {s_q:.4f}")
print(f"s_gamma = 1 - 0.3 * ({B_prime:.3f}/{L_prime:.3f}) = {s_gamma:.4f}")
print(f"s_c = ({s_q:.4f} * {N_q:.3f} - 1) / ({N_q:.3f} - 1) = {s_c:.4f}")
print("-" * 20)

# Step 7: Calculate Load Inclination Factors
m = (2 + B_prime_over_L_prime) / (1 + B_prime_over_L_prime)
# Per instruction, omitting the A'c'cot(phi) term
i_q = (1 - H_d / V_d)**m
i_gamma = (1 - H_d / V_d)**(m + 1)
i_c = (i_q * N_q - 1) / (N_q - 1)

print("--- Load Inclination Factors ---")
print(f"Exponent m = {m:.4f}")
print(f"i_q = (1 - {H_d:.2f}/{V_d:.2f})^{m:.4f} = {i_q:.4f}")
print(f"i_gamma = (1 - {H_d:.2f}/{V_d:.2f})^{m+1:.4f} = {i_gamma:.4f}")
print(f"i_c = ({i_q:.4f} * {N_q:.3f} - 1) / ({N_q:.3f} - 1) = {i_c:.4f}")
print("-" * 20)


# Step 8: Calculate Design Bearing Resistance
# The general formula is: R/A' = c'*Nc*sc*ic + q'*Nq*sq*iq + 0.5*gamma'*B'*Ngamma*sgamma*igamma
term1 = c_d * N_c * s_c * i_c
term2 = q_prime * N_q * s_q * i_q
term3 = 0.5 * gamma_d * B_prime * N_gamma * s_gamma * i_gamma
q_ult = term1 + term2 + term3
R_d = q_ult * A_prime

print("--- Final Calculation of Design Resistance, R_d ---")
print("R_d/A' = (c'_d * N_c * s_c * i_c) + (q' * N_q * s_q * i_q) + (0.5 * gamma_d * B' * N_gamma * s_gamma * i_gamma)\n")

print("Cohesion term:")
print(f"= {c_d:.2f} kPa * {N_c:.3f} * {s_c:.4f} * {i_c:.4f} = {term1:.2f} kPa\n")

print("Surcharge term:")
print(f"= {q_prime:.2f} kPa * {N_q:.3f} * {s_q:.4f} * {i_q:.4f} = {term2:.2f} kPa\n")

print("Self-weight term:")
print(f"= 0.5 * {gamma_d:.2f} kN/m^3 * {B_prime:.4f} m * {N_gamma:.3f} * {s_gamma:.4f} * {i_gamma:.4f} = {term3:.2f} kPa\n")

print(f"Design bearing capacity q_ult = {term1:.2f} + {term2:.2f} + {term3:.2f} = {q_ult:.2f} kPa\n")

print(f"Design resistance R_d = q_ult * A' = {q_ult:.2f} kPa * {A_prime:.4f} m^2 = {R_d:.1f} kN")

print("\nThe final answer for the design resistance is:")
print(f"{R_d:.1f} kN")
<<<5336.1>>>