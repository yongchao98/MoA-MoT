import math

# --- 1. Initial Parameters and Factors ---
# Characteristic Loads (kN)
Gk = 1000.0  # Permanent vertical
Qk_v = 1500.0 # Variable vertical
Qk_h = 300.0  # Variable horizontal

# Geometry (m)
h_col = 2.0      # Height of load application above footing
d_footing = 0.75 # Thickness of footing / depth of foundation base
B = 2.3          # Footing width (B=L), determined by equating Vd/A' = q_Rd
print(f"The required footing dimension is found to be B = L = {B} m.\n")

# Material Properties
gamma_concrete = 24.0 # Unit weight of concrete (kN/m^3)
phi_k_deg = 35.0      # Characteristic angle of shearing resistance (degrees)
c_prime_k = 0.0       # Characteristic effective cohesion (kPa)
gamma_soil = 20.0     # Unit weight of soil (kN/m^3)

# Partial Factors for ULS (DA1-C1)
gamma_G = 1.35  # Permanent actions
gamma_Q = 1.50  # Variable actions
gamma_phi = 1.0 # Angle of shearing resistance
gamma_c = 1.0   # Effective cohesion
gamma_gamma = 1.0 # Soil weight density

# --- 2. Calculate Design Actions ---
# Self-weight of footing base
Wk_footing = (B * B * d_footing) * gamma_concrete
Gk_total = Gk + Wk_footing

# Design loads
Vd = gamma_G * Gk_total + gamma_Q * Qk_v
Hd = gamma_Q * Qk_h
Md = Hd * (h_col + d_footing)

print("--- Design Actions and Effects ---")
print(f"Total characteristic permanent load (Gk_total) = {Gk:.1f} + {Wk_footing:.2f} = {Gk_total:.2f} kN")
print(f"Design vertical load (Vd) = {gamma_G} * {Gk_total:.2f} + {gamma_Q} * {Qk_v:.1f} = {Vd:.2f} kN")
print(f"Design horizontal load (Hd) = {gamma_Q} * {Qk_h:.1f} = {Hd:.2f} kN")
print(f"Design moment at base (Md) = {Hd:.2f} * ({h_col} + {d_footing}) = {Md:.2f} kNm")

# --- 3. Effective Footing Dimensions ---
e = Md / Vd
B_prime = B - 2 * e
L_prime = B # No eccentricity in L-direction
A_prime = B_prime * L_prime

print("\n--- Effective Footing Dimensions ---")
print(f"Eccentricity (e) = {Md:.2f} / {Vd:.2f} = {e:.4f} m")
print(f"Effective width (B') = {B:.1f} - 2 * {e:.4f} = {B_prime:.4f} m")
print(f"Effective length (L') = {L_prime:.1f} m")
print(f"Effective area (A') = {B_prime:.4f} * {L_prime:.1f} = {A_prime:.4f} m^2")

# --- 4. Bearing Resistance Calculation ---
# Design soil parameters
phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi))
c_prime_d = c_prime_k / gamma_c
gamma_d = gamma_soil / gamma_gamma
phi_d_rad = math.radians(phi_d_deg)

# Bearing capacity factors (EN 1997-1 Annex D)
Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45 + phi_d_deg / 2)))**2
Nc = (Nq - 1) / math.tan(phi_d_rad)
Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)

# Overburden pressure
q_prime = gamma_d * d_footing

# Shape factors
sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
sc = (sq * Nq - 1) / (Nq - 1)
sgamma = 1 - 0.3 * (B_prime / L_prime)

# Load inclination factors
# Exponent m
m = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
# Factors
iq = (1 - Hd / Vd)**m
ic = iq - (1 - iq) / (Nc * math.tan(phi_d_rad))
igamma = (1 - Hd / Vd)**(m + 1)

# Base and ground inclination factors are 1.0 as base and ground are horizontal.
bc = bq = bgamma = 1.0
gc = gq = ggamma = 1.0

# Bearing resistance formula terms
# q_Rd = (c'*Nc*sc*bc*ic) + (q'*Nq*sq*bq*iq) + (0.5*gamma'*B'*Ngamma*sgamma*bgamma*igamma)
term_c = c_prime_d * Nc * sc * bc * ic
term_q = q_prime * Nq * sq * bq * iq
term_gamma = 0.5 * gamma_d * B_prime * Ngamma * sgamma * bgamma * igamma
q_Rd = term_c + term_q + term_gamma

print("\n--- ULS Design Bearing Resistance Calculation (q_Rd) ---")
print("q_Rd = c'*Nc*sc*ic + q'*Nq*sq*iq + 0.5*gamma'*B'*N_gamma*s_gamma*i_gamma\n")

print("Parameters for the equation:")
print(f"  c' = {c_prime_d:.2f} kPa")
print(f"  q' = {gamma_d:.1f} * {d_footing} = {q_prime:.2f} kPa")
print(f"  gamma' = {gamma_d:.1f} kN/m^3")
print(f"  B' = {B_prime:.4f} m")
print(f"  phi'd = {phi_d_deg:.1f} degrees")
print(f"  Nc = {Nc:.3f}")
print(f"  Nq = {Nq:.3f}")
print(f"  N_gamma = {Ngamma:.3f}")
print(f"  sc = {sc:.4f}")
print(f"  sq = {sq:.4f}")
print(f"  s_gamma = {sgamma:.4f}")
print(f"  ic = {ic:.4f}")
print(f"  iq = {iq:.4f}")
print(f"  i_gamma = {igamma:.4f}\n")

print("Calculating the terms:")
print(f"Cohesion term = {c_prime_d:.2f} * {Nc:.3f} * {sc:.4f} * {ic:.4f} = {term_c:.2f} kPa")
print(f"Surcharge term = {q_prime:.2f} * {Nq:.3f} * {sq:.4f} * {iq:.4f} = {term_q:.2f} kPa")
print(f"Self-weight term = 0.5 * {gamma_d:.1f} * {B_prime:.4f} * {Ngamma:.3f} * {sgamma:.4f} * {igamma:.4f} = {term_gamma:.2f} kPa\n")
print(f"Final ULS Design Bearing Resistance (q_Rd):")
print(f"q_Rd = {term_c:.2f} + {term_q:.2f} + {term_gamma:.2f} = {q_Rd:.2f} kN/m^2")
