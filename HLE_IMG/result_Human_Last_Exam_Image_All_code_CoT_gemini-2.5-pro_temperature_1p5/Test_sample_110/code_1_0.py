import math

# Step 1: Define Material and Section Properties
# --- Given Data ---
# Geometry
L = 4000.0  # Storey height, mm
h_col = 400.0 # Column height, mm
b_col = 400.0 # Column width, mm
cover = 30.0 # Concrete cover, mm

# Concrete: C40/50
f_ck = 40.0  # Characteristic cylinder strength, N/mm^2

# Structural Steel: UC 254x254x132 S355
f_y = 345.0  # Yield strength, N/mm^2
A_a = 168.0 * 100 # Area, mm^2 (168 cm^2)
I_ay = 22500.0 * 1e4 # Second moment of area about y-y axis, mm^4 (22500 cm^4)
I_az = 7530.0 * 1e4 # Second moment of area about z-z axis, mm^4 (7530 cm^4)

# Reinforcement Steel: 4xT16 S500
num_rebars = 4
d_rebar = 16.0 # Rebar diameter, mm
f_sk = 500.0 # Characteristic yield strength, N/mm^2

# --- Eurocode Parameters ---
# Partial safety factors
gamma_a = 1.0 # for structural steel
gamma_c = 1.5 # for concrete
gamma_s = 1.15 # for reinforcing steel

# Modulus of Elasticity
E_a = 210000.0 # N/mm^2 for structural steel
E_s = 200000.0 # N/mm^2 for reinforcing steel

# Stiffness calculation coefficients
K_o = 0.9 # Calibration factor
K_e = 0.6 # Concrete contribution factor
alpha_cc = 0.85 # Factor for long-term effects on concrete strength

# Step 2: Calculate Component Areas and Strengths
# Design strengths
f_yd = f_y / gamma_a
f_cd = alpha_cc * f_ck / gamma_c
f_sd = f_sk / gamma_s

# Area of reinforcing steel
A_s_one_bar = math.pi * (d_rebar / 2)**2
A_s = num_rebars * A_s_one_bar

# Area of concrete
A_total = h_col * b_col
A_c = A_total - A_a - A_s

# Secant modulus of elasticity for concrete
E_cm = 22 * ((f_ck + 8) / 10)**0.3

# Step 3: Determine Governing Axis for Buckling
# Second moment of area for concrete section
I_c = (b_col * h_col**3) / 12

# Second moment of area for rebars (symmetric placement)
# Distance from column center to rebar center
dist_rebar_center = h_col / 2 - cover - d_rebar / 2
# Inertia of rebars using parallel axis theorem (self-inertia is negligible)
I_s = num_rebars * A_s_one_bar * dist_rebar_center**2

# Effective flexural stiffness (EI_eff)
# About y-y axis
EI_eff_y = K_o * (E_a * I_ay + E_s * I_s + K_e * E_cm * I_c)
# About z-z axis
EI_eff_z = K_o * (E_a * I_az + E_s * I_s + K_e * E_cm * I_c)

# The governing stiffness is the smaller one
EI_eff = min(EI_eff_y, EI_eff_z)

# Step 4: Calculate Plastic and Critical Loads
# Plastic resistance of cross-section (Design)
N_pl_Rd = (A_a * f_yd) + (A_c * f_cd) + (A_s * f_sd)

# Plastic resistance of cross-section (Characteristic) for slenderness calculation
N_pl_Rk = (A_a * f_y) + (alpha_cc * A_c * f_ck) + (A_s * f_sk)

# Effective length (assuming pinned-pinned column)
L_e = 1.0 * L

# Elastic critical buckling load
N_cr = (math.pi**2 * EI_eff) / L_e**2

# Step 5: Determine Buckling Reduction Factor
# Relative slenderness
lambda_bar = math.sqrt(N_pl_Rk / N_cr)

# Buckling curve 'c' for encased I-section (weak axis)
alpha_buckling = 0.49 # Imperfection factor for curve 'c'

# Intermediate value Phi for calculating reduction factor
Phi = 0.5 * (1 + alpha_buckling * (lambda_bar - 0.2) + lambda_bar**2)

# Buckling reduction factor
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
# Ensure chi is not greater than 1
chi = min(chi, 1.0)

# Step 6: Calculate Final Buckling Resistance
# Final buckling resistance
N_b_Rd = chi * N_pl_Rd

# --- Output the results ---
print("--- Calculation of Buckling Resistance ---\n")
print(f"Plastic Resistance (N_pl,Rd): {N_pl_Rd/1000:.1f} kN")
print(f"Buckling Reduction Factor (chi): {chi:.4f}")
print("\n--- Final Buckling Resistance Calculation ---")
print(f"N_b,Rd = chi * N_pl,Rd")
print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd/1000:.1f} kN")

final_answer_kN = N_b_Rd / 1000
print(f"\nThe buckling resistance of the column is: {final_answer_kN:.1f} kN")
print(f"\n<<< {final_answer_kN:.1f} >>>")