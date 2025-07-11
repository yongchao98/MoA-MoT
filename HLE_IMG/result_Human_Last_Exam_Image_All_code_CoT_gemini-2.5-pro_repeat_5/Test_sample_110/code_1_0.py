import math

# --- 1. Define Material and Section Properties ---
# Storey height (m)
L = 4.0
# Effective length factor (assuming pinned ends)
k = 1.0
L_cr = L * k * 1000  # in mm

# Structural Steel (S355 - UC 254x254x132)
f_y = 345.0  # Yield strength (N/mm^2)
A_a = 168.0 * 100  # Area (mm^2)
I_ay = 22500.0 * 1e4  # Second moment of area about y-y (mm^4)
I_az = 7530.0 * 1e4  # Second moment of area about z-z (mm^4)
E_a = 210000.0  # Modulus of Elasticity (N/mm^2)
gamma_a = 1.0  # Partial safety factor

# Concrete (C40/50)
f_ck = 40.0  # Characteristic cylinder strength (N/mm^2)
h_c = 400.0  # Concrete section height (mm)
b_c = 400.0  # Concrete section width (mm)
E_cm = 35000.0 # Secant modulus of elasticity (N/mm^2)
gamma_c = 1.5  # Partial safety factor
alpha_c = 0.85 # Factor for concrete strength in compression

# Reinforcing Steel (S500 - 4xT16)
f_sk = 500.0  # Characteristic yield strength (N/mm^2)
d_s = 16.0  # Rebar diameter (mm)
num_s = 4  # Number of rebars
E_s = 200000.0  # Modulus of Elasticity (N/mm^2)
gamma_s = 1.15  # Partial safety factor
cover = 30.0  # Concrete cover (mm)

print("--- Step 1: Input Properties ---")
print(f"Column effective length, L_cr = {L_cr:.0f} mm")
print(f"Steel yield strength, f_y = {f_y} MPa")
print(f"Concrete characteristic strength, f_ck = {f_ck} MPa")
print(f"Rebar characteristic strength, f_sk = {f_sk} MPa")
print("-" * 30 + "\n")

# --- 2. Calculate Cross-Sectional Properties ---
# Rebar properties
A_s1 = math.pi * (d_s / 2)**2  # Area of one rebar
A_s = num_s * A_s1
# Distance from global centroid to rebar centroid
d_rebar = h_c / 2 - cover - d_s / 2

# Moments of inertia for reinforcement (I_s_centroidal is negligible)
# Using Parallel Axis Theorem: I = I_centroidal + A*d^2
I_sy = num_s * (A_s1 * d_rebar**2)
I_sz = num_s * (A_s1 * d_rebar**2)

# Concrete properties
A_c_gross = h_c * b_c
A_c = A_c_gross - A_a - A_s
I_c_gross = h_c**4 / 12
I_cy = I_c_gross - I_ay - I_sy
I_cz = I_c_gross - I_az - I_sz

print("--- Step 2: Sectional Properties ---")
print(f"Total rebar area, A_s = {A_s:.2f} mm^2")
print(f"Net concrete area, A_c = {A_c:.2f} mm^2")
print(f"Rebar moment of inertia about z-axis, I_sz = {I_sz:.2e} mm^4")
print(f"Concrete moment of inertia about z-axis, I_cz = {I_cz:.2e} mm^4")
print(f"Steel section moment of inertia about z-axis, I_az = {I_az:.2e} mm^4")
print("-" * 30 + "\n")


# --- 3. Calculate Plastic Resistance (N_pl,Rd) ---
f_yd = f_y / gamma_a
f_cd = f_ck / gamma_c
f_sd = f_sk / gamma_s
N_pl_Rd_a = A_a * f_yd
N_pl_Rd_c = alpha_c * A_c * f_cd
N_pl_Rd_s = A_s * f_sd
N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s
N_pl_Rd_kN = N_pl_Rd / 1000

print("--- Step 3: Plastic Resistance (N_pl,Rd) ---")
print(f"Design Plastic Resistance, N_pl,Rd = {N_pl_Rd_kN:.1f} kN")
print("-" * 30 + "\n")


# --- 4. Determine Effective Flexural Stiffness ((EI)_eff) ---
# Check weaker axis (lower I value for steel section is I_az)
K_e = 0.6  # Correction factor for E_cm
EI_eff_y = K_e * E_a * I_ay + E_s * I_sy + K_e * E_cm * I_cy
EI_eff_z = K_e * E_a * I_az + E_s * I_sz + K_e * E_cm * I_cz
EI_eff = min(EI_eff_y, EI_eff_z)

print("--- Step 4: Effective Flexural Stiffness ((EI)_eff) ---")
print(f"Effective stiffness about y-axis, (EI)_eff,y = {EI_eff_y:.2e} Nmm^2")
print(f"Effective stiffness about z-axis, (EI)_eff,z = {EI_eff_z:.2e} Nmm^2")
print(f"Governing (minimum) stiffness, (EI)_eff = {EI_eff:.2e} Nmm^2")
print("-" * 30 + "\n")


# --- 5. Calculate Critical Buckling Load (N_cr) ---
N_cr = (math.pi**2 * EI_eff) / (L_cr**2)
N_cr_kN = N_cr / 1000

print("--- Step 5: Critical Buckling Load (N_cr) ---")
print(f"Critical Buckling Load, N_cr = {N_cr_kN:.1f} kN")
print("-" * 30 + "\n")


# --- 6. Calculate Relative Slenderness (lambda_bar) ---
# N_pl,Rk uses characteristic strengths. A factor of 0.85 is commonly applied
# to the concrete contribution for encased sections, consistent with design practice.
N_pl_Rk = A_a * f_y + alpha_c * A_c * f_ck + A_s * f_sk
lambda_bar = math.sqrt(N_pl_Rk / N_cr)

print("--- Step 6: Relative Slenderness (lambda_bar) ---")
print(f"Characteristic Plastic Resistance, N_pl,Rk = {N_pl_Rk/1000:.1f} kN")
print(f"Relative Slenderness, lambda_bar = {lambda_bar:.3f}")
print("-" * 30 + "\n")


# --- 7. Determine Reduction Factor (chi) ---
# For encased I-section, weak axis buckling -> Buckling Curve 'b'
alpha = 0.34  # Imperfection factor for buckling curve 'b'
Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
if chi > 1.0:
    chi = 1.0

print("--- Step 7: Reduction Factor (chi) ---")
print(f"Using Buckling Curve 'b', imperfection factor alpha = {alpha}")
print(f"Intermediate value, Phi = {Phi:.3f}")
print(f"Buckling Reduction Factor, chi = {chi:.3f}")
print("-" * 30 + "\n")


# --- 8. Calculate Buckling Resistance (N_b,Rd) ---
N_b_Rd = chi * N_pl_Rd
N_b_Rd_kN = N_b_Rd / 1000

print("--- Step 8: Final Buckling Resistance (N_b,Rd) ---")
print(f"The final design equation is N_b,Rd = chi * N_pl,Rd")
print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd_kN:.1f} kN")
print(f"Design Buckling Resistance, N_b,Rd = {N_b_Rd_kN:.1f} kN")
print("-" * 30 + "\n")

print(f"The buckling resistance of the column is {N_b_Rd_kN:.1f} kN.")
<<<7917.4>>>