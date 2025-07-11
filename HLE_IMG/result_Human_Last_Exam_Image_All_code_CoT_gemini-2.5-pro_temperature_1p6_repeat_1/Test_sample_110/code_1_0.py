import math

# Step 1: Define Material and Geometric Properties
# --- Steel Section (UC 254 x 254 x 132, S355) ---
f_y = 345.0  # N/mm^2, Yield strength
A_a = 168.0 * 100  # mm^2, Area (168 cm^2)
I_ay = 22500.0 * 1e4 # mm^4, Second moment of area about y-y (22500 cm^4)
I_az = 7530.0 * 1e4  # mm^4, Second moment of area about z-z (7530 cm^4)
E_a = 210000.0 # N/mm^2, Modulus of Elasticity for steel

# --- Concrete (C40/50) ---
f_ck = 40.0  # N/mm^2, Characteristic cylinder strength
# E_cm from EN 1992-1-1 Table 3.1: E_cm = 22 * ((f_ck + 8)/10)^0.3
E_cm = 22 * ((f_ck + 8) / 10)**0.3

# --- Reinforcement (4 no. T16, S500) ---
n_s = 4  # Number of rebars
d_s = 16.0  # mm, Diameter of rebar
f_sk = 500.0  # N/mm^2, Characteristic yield strength of rebar
E_s = 200000.0 # N/mm^2, Modulus of Elasticity for rebar

# --- Column Geometry ---
h_c = 400.0  # mm, Column height
b_c = 400.0  # mm, Column width
cover = 30.0  # mm
L = 4000.0 # mm, Column storey height
L_cr = 1.0 * L # Effective length for pinned-pinned ends

# --- Partial Safety Factors (from Eurocodes) ---
gamma_c = 1.5
gamma_s = 1.15
gamma_a = 1.0 # For structural steel in compression
gamma_M1 = 1.0 # For buckling resistance

# --- Derived design strengths ---
f_cd = f_ck / gamma_c
f_sd = f_sk / gamma_s
f_yd = f_y / gamma_a

print("--- Step 1: Material and Geometric Properties ---")
print(f"Concrete C40/50: f_ck = {f_ck} MPa, E_cm = {E_cm:.0f} MPa")
print(f"Steel UC254x254x132 S355: f_y = {f_y} MPa, A_a = {A_a} mm^2, I_ay = {I_ay:.0f} mm^4, I_az = {I_az:.0f} mm^4")
print(f"Rebar S500: {n_s} no. T{int(d_s)}, f_sk = {f_sk} MPa")
print(f"Column: {h_c}x{b_c} mm, Length L_cr = {L_cr} mm")
print("-" * 50)

# Step 2: Calculate Section Areas
A_s_one = math.pi * (d_s / 2)**2
A_s = n_s * A_s_one
A_gross_c = h_c * b_c
# Net concrete area
A_c = A_gross_c - A_a - A_s

print("--- Step 2: Sectional Areas ---")
print(f"Total reinforcement area, A_s = {A_s:.2f} mm^2")
print(f"Net concrete area, A_c = {A_c:.2f} mm^2")
print("-" * 50)

# Step 3: Calculate Plastic Design Resistance (N_pl,Rd)
# According to EN 1994-1-1 6.7.3.2(2), use alpha_c = 1.0 for encased sections
alpha_c_plastic = 1.0
N_pl_Rd = (A_a * f_yd) + (alpha_c_plastic * A_c * f_cd) + (A_s * f_sd)

# Calculate Characteristic Plastic Resistance (N_pl,Rk) for slenderness calculation
N_pl_Rk = (A_a * f_y) + (alpha_c_plastic * A_c * f_ck) + (A_s * f_sk)

print("--- Step 3: Plastic Resistance of Section ---")
print(f"Design Plastic Resistance, N_pl,Rd = {N_pl_Rd / 1000:.1f} kN")
print(f"Characteristic Plastic Resistance, N_pl,Rk = {N_pl_Rk / 1000:.1f} kN")
print("-" * 50)


# Step 4: Calculate Effective Flexural Stiffness ((EI)_eff)
# We check the weaker axis (lower I value). For the UC section, I_az < I_ay.
# The overall section is square, so the weaker axis is governed by the steel section's weaker axis.
# Buckling will be about the z-z axis.

# Second moment of area of gross concrete section
I_cz = (b_c * h_c**3) / 12

# Second moment of area of reinforcement about z-z axis
dist_rebar_z = h_c / 2 - cover - d_s / 2
I_sz = A_s * dist_rebar_z**2 # Approximation ignoring self-I of rebars

# Per problem, ignore long-term effects, so creep coefficient is 0.
# The factor 0.6 in the standard EC4 formula for (EI)eff is for creep, so we take it as 1.0.
K_e = 1.0 
K_o = 0.9 # Calibration factor from EC4 simplified method

EI_eff_z = K_o * (E_a * I_az + E_s * I_sz + K_e * E_cm * I_cz)

print("--- Step 4: Effective Flexural Stiffness ((EI)_eff) ---")
print(f"Buckling is critical about the z-z axis.")
print(f"Stiffness of steel section (E_a*I_az): {E_a*I_az:.2e} Nmm^2")
print(f"Stiffness of concrete (K_e*E_cm*I_cz): {K_e*E_cm*I_cz:.2e} Nmm^2")
print(f"Stiffness of rebar (E_s*I_sz): {E_s*I_sz:.2e} Nmm^2")
print(f"Total Effective Flexural Stiffness, (EI)_eff = {EI_eff_z:.2e} Nmm^2")
print("-" * 50)


# Step 5: Calculate Critical Buckling Load (N_cr)
N_cr = (math.pi**2 * EI_eff_z) / (L_cr**2)

print("--- Step 5: Critical Buckling Load (N_cr) ---")
print(f"N_cr = (pi^2 * (EI)_eff) / L_cr^2")
print(f"N_cr = {N_cr / 1000:.1f} kN")
print("-" * 50)

# Step 6: Calculate Relative Slenderness (lambda_bar)
lambda_bar = math.sqrt(N_pl_Rk / N_cr)

print("--- Step 6: Relative Slenderness (lambda_bar) ---")
print(f"lambda_bar = sqrt(N_pl,Rk / N_cr)")
print(f"lambda_bar = sqrt({N_pl_Rk/1000:.1f} / {N_cr/1000:.1f}) = {lambda_bar:.4f}")
print("-" * 50)

# Step 7: Determine Buckling Reduction Factor (chi)
# From EN 1994-1-1 Table 6.5, for concrete encased I-section, buckling about z-z axis
# and concrete h/b = 400/400 = 1.0 <= 1.2, use buckling curve 'c'.
# For buckling curve 'c', imperfection factor alpha = 0.49 (from EN 1993-1-1 Table 6.1)
alpha = 0.49
print("--- Step 7: Buckling Reduction Factor (chi) ---")
print("From EN 1994-1-1 Table 6.5, for weak axis (z-z) buckling, curve 'c' applies.")
print(f"Imperfection factor for curve 'c', alpha = {alpha}")

Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
# Ensure chi is not > 1.0
chi = min(chi, 1.0)

print(f"Phi = 0.5 * [1 + {alpha} * ({lambda_bar:.4f} - 0.2) + {lambda_bar:.4f}^2] = {Phi:.4f}")
print(f"chi = 1 / (Phi + sqrt(Phi^2 - lambda_bar^2)) = {chi:.4f}")
print("-" * 50)

# Step 8: Calculate Design Buckling Resistance (N_b,Rd)
N_b_Rd = chi * N_pl_Rd

print("--- Step 8: Final Buckling Resistance (N_b,Rd) ---")
print("The design buckling resistance is calculated as:")
print("N_b,Rd = chi * N_pl,Rd")
print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd / 1000:.1f} kN")
final_answer_kN = N_b_Rd / 1000
print(f"N_b,Rd = {final_answer_kN:.1f} kN")
print("-" * 50)

print("The buckling resistance of the column is:")
print(f"{final_answer_kN:.1f} kN")
