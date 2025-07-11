import math

# Step 1: Define the given data and constants
G = 500.0  # Permanent vertical load (kN)
Q_v = 1500.0 # Variable vertical load (kN)
Q_h = 120.0  # Variable horizontal load (kN)
zeta = 0.7   # Combination factor for independent variable loads
B = 2.39     # Foundation width (m)
L = 2.39     # Foundation length (m)
z = 1.0      # Foundation depth (m)
c_k = 10.0   # Characteristic effective cohesion (kPa or kN/m^2)
phi_k_deg = 28.0 # Characteristic effective friction angle (degrees)
gamma_soil = 20.0 # Soil unit weight (kN/m^3)
# Eurocode 7 partial factors for Load Combination 1 (DA1-1)
gamma_G = 1.35 # Partial factor for permanent actions
gamma_Q = 1.5  # Partial factor for variable actions
# Partial factors for soil parameters are 1.0 for DA1-1
gamma_c = 1.0
gamma_phi = 1.0

print("This script calculates the design bearing resistance of the square pad foundation.")
print("-" * 30)
print("1. Initial Data:")
print(f"  Permanent vertical load, G = {G} kN")
print(f"  Variable vertical load, Q_v = {Q_v} kN")
print(f"  Variable horizontal load, Q_h = {Q_h} kN")
print(f"  Foundation width, B = L = {B} m")
print(f"  Foundation depth, z = {z} m")
print(f"  Effective cohesion, c'_k = {c_k} kPa")
print(f"  Effective friction angle, phi'_k = {phi_k_deg} degrees")
print(f"  Soil unit weight, gamma_soil = {gamma_soil} kN/m^3")
print("-" * 30)

# Step 2: Calculate design actions and eccentricity
V_d = gamma_G * G + gamma_Q * Q_v
H_d = gamma_Q * zeta * Q_h
M_d = H_d * z # Moment at the foundation base
e_B = M_d / V_d # Eccentricity

print("2. Design Actions (Loads):")
print(f"  Design vertical load, V_d = {gamma_G}*{G} + {gamma_Q}*{Q_v} = {V_d:.2f} kN")
print(f"  Design horizontal load, H_d = {gamma_Q}*{zeta}*{Q_h} = {H_d:.2f} kN")
print(f"  Moment at base, M_d = {H_d:.2f} * {z:.1f} = {M_d:.2f} kNm")
print(f"  Eccentricity, e = M_d / V_d = {M_d:.2f} / {V_d:.2f} = {e_B:.4f} m")
print("-" * 30)

# Step 3: Calculate effective foundation dimensions
B_prime = B - 2 * e_B
L_prime = L
A_prime = B_prime * L_prime

print("3. Effective Foundation Dimensions:")
print(f"  Effective width, B' = B - 2*e = {B:.2f} - 2*{e_B:.4f} = {B_prime:.4f} m")
print(f"  Effective length, L' = {L_prime:.2f} m")
print(f"  Effective area, A' = B' * L' = {B_prime:.4f} * {L_prime:.2f} = {A_prime:.4f} m^2")
print("-" * 30)

# Step 4: Calculate design soil parameters
c_d = c_k / gamma_c
phi_d_deg = phi_k_deg
phi_d_rad = math.radians(phi_d_deg)

# Step 5: Calculate bearing capacity factors
tan_phi = math.tan(phi_d_rad)
Nq = math.exp(math.pi * tan_phi) * (math.tan(math.radians(45 + phi_d_deg / 2)))**2
Nc = (Nq - 1) / tan_phi
Ny = 2 * (Nq - 1) * tan_phi

print("4. Bearing Capacity Factors (for phi' = 28 degrees):")
print(f"  Nq = {Nq:.2f}")
print(f"  Nc = {Nc:.2f}")
print(f"  Ny = {Ny:.2f}")
print("-" * 30)

# Step 6: Calculate shape factors
B_prime_over_L_prime = B_prime / L_prime
s_q = 1 + B_prime_over_L_prime * math.sin(phi_d_rad)
s_y = 1 - 0.3 * B_prime_over_L_prime
s_c = (s_q * Nq - 1) / (Nq - 1)

print("5. Shape Factors:")
print(f"  sq = 1 + (B'/L')*sin(phi') = 1 + {B_prime_over_L_prime:.3f}*sin({phi_d_deg}) = {s_q:.3f}")
print(f"  sy = 1 - 0.3*(B'/L') = 1 - 0.3*{B_prime_over_L_prime:.3f} = {s_y:.3f}")
print(f"  sc = (sq*Nq - 1)/(Nq - 1) = ({s_q:.3f}*{Nq:.2f} - 1)/({Nq:.2f} - 1) = {s_c:.3f}")
print("-" * 30)

# Step 7: Calculate inclination factors
m = (2 + B_prime_over_L_prime) / (1 + B_prime_over_L_prime)
# The term A'c'cot(phi) is omitted for conservativeness as stated in the problem
ratio_H_V = H_d / V_d
i_q = (1 - ratio_H_V)**m
i_y = (1 - ratio_H_V)**(m + 1)
i_c = i_q - (1 - i_q) / (Nc * tan_phi)

print("6. Inclination Factors:")
print(f"  m = (2 + B'/L')/(1 + B'/L') = {m:.3f}")
print(f"  iq = (1 - H/V)^m = (1 - {H_d:.2f}/{V_d:.2f})^{m:.3f} = {i_q:.3f}")
print(f"  iy = (1 - H/V)^(m+1) = (1 - {H_d:.2f}/{V_d:.2f})^{m+1:.3f} = {i_y:.3f}")
print(f"  ic = iq - (1 - iq)/(Nc*tan(phi')) = {i_q:.3f} - (1-{i_q:.3f})/({Nc:.2f}*tan({phi_d_deg})) = {i_c:.3f}")
print("-" * 30)

# Step 8: Calculate overburden pressure
q_prime = gamma_soil * z

print(f"7. Overburden Pressure at base, q' = {gamma_soil} * {z} = {q_prime} kPa")
print("-" * 30)

# Step 9: Calculate design resistance
# R_d = A' * (c'd * Nc * sc * ic + q' * Nq * sq * iq + 0.5 * gamma' * B' * Ny * sy * iy)
term1 = c_d * Nc * s_c * i_c
term2 = q_prime * Nq * s_q * i_q
term3 = 0.5 * gamma_soil * B_prime * Ny * s_y * i_y
q_d = term1 + term2 + term3
R_d = q_d * A_prime

print("8. Final Calculation of Design Resistance, R_d:")
print("The design resistance is calculated using the formula:")
print("R_d = A' * (c'_d*Nc*sc*ic + q'*Nq*sq*iq + 0.5*gamma'*B'*Ny*sy*iy)\n")

print("Substituting the values:")
print(f"R_d = {A_prime:.4f} * ({c_d:.1f}*{Nc:.2f}*{s_c:.3f}*{i_c:.3f} + "
      f"{q_prime:.1f}*{Nq:.2f}*{s_q:.3f}*{i_q:.3f} + "
      f"0.5*{gamma_soil:.1f}*{B_prime:.4f}*{Ny:.2f}*{s_y:.3f}*{i_y:.3f})")

print(f"\nCalculating each term of the bearing pressure:")
print(f"Cohesion term = {term1:.2f} kPa")
print(f"Surcharge term = {term2:.2f} kPa")
print(f"Self-weight term = {term3:.2f} kPa")

print(f"\nTotal design bearing pressure, q_d = {term1:.2f} + {term2:.2f} + {term3:.2f} = {q_d:.2f} kPa")
print(f"Total design resistance, R_d = q_d * A' = {q_d:.2f} * {A_prime:.4f} = {R_d:.2f} kN")
print("\nFinal Answer:")
print(f"The design resistance under load combination 1 is {R_d:.1f} kN.")