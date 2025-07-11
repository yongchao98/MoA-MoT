import math

# Step 1: Define Material and Geometric Properties
# --- Given Data ---
# Concrete C40/50
f_ck = 40.0  # N/mm^2, characteristic cylinder strength
# Steel Section UC 254x254x132 S355
f_y = 345.0  # N/mm^2, yield strength
A_a = 168.0 * 100  # mm^2, Area from cm^2
I_az = 7530.0 * 10000  # mm^4, Second moment of area about z-z from cm^4
# Reinforcing Steel T16, S500
f_sk = 500.0  # N/mm^2, characteristic yield strength
num_rebars = 4
d_rebar = 16.0  # mm
# Column Geometry
H_col = 400.0  # mm
B_col = 400.0  # mm
cover = 30.0  # mm
storey_height = 4000.0  # mm
# Assume pinned-pinned connection, effective length factor k=1.0
L_cr = 1.0 * storey_height # mm, effective buckling length

# --- Standard Properties & Factors (Eurocode) ---
E_a = 210000.0  # N/mm^2, Modulus of elasticity for structural steel
E_s = 200000.0  # N/mm^2, Modulus of elasticity for reinforcing steel
E_cm = 35000.0  # N/mm^2, Secant modulus of elasticity of C40/50 concrete (from EC2 Table 3.1)
gamma_a = 1.0  # Partial factor for structural steel
gamma_c = 1.5  # Partial factor for concrete
gamma_s = 1.15 # Partial factor for reinforcing steel
K_e_c = 0.6    # Coefficient for concrete contribution to stiffness (short-term)
alpha = 0.49   # Imperfection factor for buckling curve 'c' (for encased I-section about z-z axis)

# Step 2: Calculate Cross-Sectional Properties
A_s_one_rebar = math.pi * (d_rebar / 2)**2
A_s = num_rebars * A_s_one_rebar  # Total area of reinforcement
A_c_gross = H_col * B_col
A_c = A_c_gross - A_a - A_s # Net area of concrete

# Second moment of area of reinforcement about z-z axis
# Distance from column centroid to rebar centroid
rebar_dist_z = (H_col / 2) - cover - (d_rebar / 2)
# Using Steiner's theorem (I = I_self + Ad^2), I_self is ignored for rebars
I_sz = A_s * (rebar_dist_z**2)
# Second moment of area of gross concrete section about z-z axis
I_cz = B_col * (H_col**3) / 12

print("--- Step-by-step Calculation for Buckling Resistance ---\n")

# Step 3: Calculate Plastic Resistance (N_pl,Rd)
print("Step 1: Calculate Design Plastic Resistance (N_pl,Rd)")
print("Formula: N_pl,Rd = A_a*f_y/γ_a + 0.85*A_c*f_ck/γ_c + A_s*f_sk/γ_s")
N_pl_Rd_a = A_a * f_y / gamma_a
N_pl_Rd_c = 0.85 * A_c * f_ck / gamma_c
N_pl_Rd_s = A_s * f_sk / gamma_s
N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s
print(f"N_pl,Rd = ({A_a:.1f}*{f_y:.1f}/{gamma_a:.1f}) + (0.85*{A_c:.1f}*{f_ck:.1f}/{gamma_c:.1f}) + ({A_s:.1f}*{f_sk:.1f}/{gamma_s:.2f})")
print(f"N_pl,Rd = {N_pl_Rd_a:.0f} + {N_pl_Rd_c:.0f} + {N_pl_Rd_s:.0f} = {N_pl_Rd:.0f} N\n")

# Step 4: Calculate Effective Flexural Stiffness ((EI)_eff,z)
print("Step 2: Calculate Effective Flexural Stiffness ((EI)_eff,z)")
print("Formula: (EI)_eff = E_a*I_az + E_s*I_sz + K_e,c*E_cm*I_cz")
EI_eff_a = E_a * I_az
EI_eff_s = E_s * I_sz
EI_eff_c = K_e_c * E_cm * I_cz
EI_eff_z = EI_eff_a + EI_eff_s + EI_eff_c
print(f"(EI)_eff,z = ({E_a:.0f}*{I_az:.0f}) + ({E_s:.0f}*{I_sz:.0f}) + ({K_e_c:.1f}*{E_cm:.0f}*{I_cz:.0f})")
print(f"(EI)_eff,z = {EI_eff_a:.2e} + {EI_eff_s:.2e} + {EI_eff_c:.2e} = {EI_eff_z:.2e} N.mm^2\n")

# Step 5: Determine Buckling Parameters
# Characteristic plastic resistance (N_pl,Rk)
print("Step 3: Calculate Characteristic Plastic Resistance (N_pl,Rk)")
print("Formula: N_pl,Rk = A_a*f_y + 0.85*A_c*f_ck + A_s*f_sk")
N_pl_Rk_a = A_a * f_y
N_pl_Rk_c = 0.85 * A_c * f_ck
N_pl_Rk_s = A_s * f_sk
N_pl_Rk = N_pl_Rk_a + N_pl_Rk_c + N_pl_Rk_s
print(f"N_pl,Rk = ({A_a:.1f}*{f_y:.1f}) + (0.85*{A_c:.1f}*{f_ck:.1f}) + ({A_s:.1f}*{f_sk:.1f})")
print(f"N_pl,Rk = {N_pl_Rk_a:.0f} + {N_pl_Rk_c:.0f} + {N_pl_Rk_s:.0f} = {N_pl_Rk:.0f} N\n")

# Elastic critical buckling load (N_cr)
print("Step 4: Calculate Elastic Critical Buckling Load (N_cr,z)")
print("Formula: N_cr,z = π² * (EI)_eff,z / L_cr²")
N_cr_z = (math.pi**2 * EI_eff_z) / (L_cr**2)
print(f"N_cr,z = π² * {EI_eff_z:.2e} / {L_cr:.0f}² = {N_cr_z:.0f} N\n")

# Relative slenderness (lambda_bar)
print("Step 5: Calculate Relative Slenderness (λ_bar)")
print("Formula: λ_bar = sqrt(N_pl,Rk / N_cr,z)")
lambda_bar = math.sqrt(N_pl_Rk / N_cr_z)
print(f"λ_bar = sqrt({N_pl_Rk:.0f} / {N_cr_z:.0f}) = {lambda_bar:.3f}\n")

# Step 6: Calculate Buckling Reduction Factor (χ)
print("Step 6: Calculate Buckling Reduction Factor (χ)")
print("Using buckling curve 'c' (α = 0.49)")
print("Formula: Φ = 0.5 * [1 + α * (λ_bar - 0.2) + λ_bar²]")
Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
print(f"Φ = 0.5 * [1 + {alpha:.2f} * ({lambda_bar:.3f} - 0.2) + {lambda_bar:.3f}²] = {Phi:.3f}")
print("Formula: χ = 1 / (Φ + sqrt(Φ² - λ_bar²))")
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
print(f"χ = 1 / ({Phi:.3f} + sqrt({Phi:.3f}² - {lambda_bar:.3f}²)) = {chi:.3f}\n")

# Step 7: Calculate Buckling Resistance (N_b,Rd)
print("Step 7: Calculate Final Buckling Resistance (N_b,Rd)")
print("Formula: N_b,Rd = χ * N_pl,Rd")
N_b_Rd = chi * N_pl_Rd
N_b_Rd_kN = N_b_Rd / 1000.0
print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd:.0f} = {N_b_Rd:.0f} N")
print(f"\nThe buckling resistance of the column is {N_b_Rd_kN:.1f} kN.")

print(f"\nFinal Answer: {N_b_Rd_kN:.1f}")