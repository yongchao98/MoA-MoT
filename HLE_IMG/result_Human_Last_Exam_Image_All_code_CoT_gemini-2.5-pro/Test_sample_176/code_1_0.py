import math

# Step 1: Define Input Parameters
G_k = 500.0  # Permanent vertical load (kN)
Q_v_k = 1500.0 # Variable vertical load (kN)
Q_h_k = 120.0  # Variable horizontal load (kN)
zeta = 0.7   # Combination factor for independent variable loads

B = 2.39     # Foundation width (m)
L = 2.39     # Foundation length (m)
d = 1.0      # Foundation depth (m)
h_f = 1.0    # Foundation thickness, assumed equal to depth (m)

c_k = 10.0     # Characteristic effective cohesion (kPa or kN/m^2)
phi_k_deg = 28.0 # Characteristic effective angle of friction (degrees)
gamma_soil = 20.0  # Soil unit weight (kN/m^3)
gamma_concrete = 24.0 # Concrete unit weight (kN/m^3)

# Partial factors for Eurocode 7, DA1-1
gamma_G = 1.35
gamma_Q = 1.5
# Material partial factors for DA1-1 are 1.0
gamma_c = 1.0
gamma_phi = 1.0
gamma_gamma = 1.0

# Step 2: Calculate Design Actions
# Foundation self-weight (permanent load)
W_f = B * L * h_f * gamma_concrete
# Design vertical load
V_d = gamma_G * (G_k + W_f) + gamma_Q * Q_v_k
# Design horizontal load
H_d = gamma_Q * zeta * Q_h_k

# Step 3: Determine Effective Foundation Dimensions
# Moment at the base
M_d = H_d * h_f
# Eccentricity
e = M_d / V_d
# Effective dimensions
B_prime = B - 2 * e
L_prime = L
A_prime = B_prime * L_prime

# Step 4: Calculate Factors
# Design soil parameters (for DA1-1, design = characteristic)
c_d = c_k / gamma_c
phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi))
phi_d_rad = math.radians(phi_d_deg)

# Overburden pressure at foundation base
q_prime = gamma_soil * d * gamma_gamma

# Bearing capacity factors
tan_phi = math.tan(phi_d_rad)
N_q = math.exp(math.pi * tan_phi) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
N_c = (N_q - 1) / tan_phi
N_gamma = 2 * (N_q - 1) * tan_phi

# Shape factors
s_q = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
s_gamma = 1 - 0.3 * (B_prime / L_prime)
s_c = (s_q * N_q - 1) / (N_q - 1)

# Inclination factors
# Exponent m
m = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
# As per problem, omit A'c'cot(phi) term
i_q = (1 - H_d / V_d)**m
i_gamma = (1 - H_d / V_d)**(m + 1)
i_c = i_q - (1 - i_q) / (N_c * tan_phi)

# Step 5: Calculate Design Bearing Resistance
# Calculate each term of the bearing capacity (in kPa)
term_c = c_d * N_c * s_c * i_c
term_q = q_prime * N_q * s_q * i_q
term_gamma = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * i_gamma

# Total ultimate bearing pressure (kPa)
q_ult = term_c + term_q + term_gamma

# Design resistance (kN)
R_d = A_prime * q_ult

# Print the final equation with all calculated values
print("Calculation of Design Bearing Resistance (Rd)")
print("----------------------------------------------")
print(f"The general equation is: Rd = A' * (c'*Nc*sc*ic + q'*Nq*sq*iq + 0.5*gamma*B'*N_gamma*s_gamma*i_gamma)\n")

print("Calculated values for each parameter:")
print(f"A' = {A_prime:.3f} m^2")
print(f"B' = {B_prime:.3f} m")
print(f"c'd = {c_d:.1f} kPa")
print(f"q' = {q_prime:.1f} kPa")
print(f"gamma = {gamma_soil:.1f} kN/m^3")
print(f"Nc = {N_c:.2f}, Nq = {N_q:.2f}, N_gamma = {N_gamma:.2f}")
print(f"sc = {s_c:.3f}, sq = {s_q:.3f}, s_gamma = {s_gamma:.3f}")
print(f"ic = {i_c:.3f}, iq = {i_q:.3f}, i_gamma = {i_gamma:.3f}\n")

print("Substituting the values into the equation:")
print(f"Rd = {A_prime:.3f} * ({c_d:.1f} * {N_c:.2f} * {s_c:.3f} * {i_c:.3f} + "
      f"{q_prime:.1f} * {N_q:.2f} * {s_q:.3f} * {i_q:.3f} + "
      f"0.5 * {gamma_soil:.1f} * {B_prime:.3f} * {N_gamma:.2f} * {s_gamma:.3f} * {i_gamma:.3f})")

print(f"Rd = {A_prime:.3f} * ({term_c:.1f} + {term_q:.1f} + {term_gamma:.1f})")
print(f"Rd = {A_prime:.3f} * {q_ult:.1f}")
print(f"Rd = {R_d:.1f} kN\n")
print("The design resistance is:")
print(f"{R_d:.1f}")