import math

# Step 1: Define Inputs & Constants (in N and mm)
print("Step 1: Input parameters and constants")
# --- Column properties
storey_height = 4000.0  # mm
# Assume pinned-pinned ends based on "typical storey"
buckling_length_factor = 1.0
L_cr = storey_height * buckling_length_factor

# --- Material Properties
f_ck = 40.0  # N/mm^2 (C40/50 concrete)
f_y = 345.0  # N/mm^2 (S355 steel)
f_sk = 500.0 # N/mm^2 (S500 rebar)
E_a = 210000.0  # N/mm^2 (Modulus of elasticity for steel)
E_s = 210000.0  # N/mm^2 (Modulus of elasticity for rebar)

# --- Sectional Properties
# Concrete section
h_c = 400.0  # mm
b_c = 400.0  # mm
# Steel section (UC 254 x 254 x 132)
A_a = 168.0 * 100  # cm^2 to mm^2
I_az = 7530.0 * 10**4 # cm^4 to mm^4 (weak axis)
# Reinforcement (4 no. T16)
num_rebars = 4
d_rebar = 16.0  # mm
cover = 30.0    # mm

# --- Eurocode Factors
gamma_c = 1.5
gamma_s = 1.15
gamma_M1 = 1.0
alpha_cc = 0.85
# Stiffness reduction factor for concrete from EN 1994-1-1, 6.7.3.3(3)
K_e = 0.6
print(f"Storey height = {storey_height} mm, Buckling length L_cr = {L_cr} mm\n")


# Step 2: Calculate Material Design Strengths
print("Step 2: Material Design Strengths")
f_cd = alpha_cc * f_ck / gamma_c
f_yd = f_y / gamma_M1
f_sd = f_sk / gamma_s
print(f"f_cd = {f_cd:.2f} N/mm^2")
print(f"f_yd = {f_yd:.2f} N/mm^2")
print(f"f_sd = {f_sd:.2f} N/mm^2\n")

# Step 3: Calculate Sectional Areas
print("Step 3: Sectional Areas")
A_c_gross = h_c * b_c
A_c = A_c_gross - A_a
A_s_one_rebar = math.pi * (d_rebar / 2)**2
A_s = num_rebars * A_s_one_rebar
print(f"Area of steel section A_a = {A_a:.2f} mm^2")
print(f"Area of concrete A_c = {A_c:.2f} mm^2")
print(f"Area of reinforcement A_s = {A_s:.2f} mm^2\n")

# Step 4: Calculate Plastic Compressive Resistance (N_pl,Rd)
print("Step 4: Plastic Compressive Resistance")
Npl_Rd_a = A_a * f_yd
Npl_Rd_c = A_c * f_cd
Npl_Rd_s = A_s * f_sd
N_pl_Rd = Npl_Rd_a + Npl_Rd_c + Npl_Rd_s
print(f"N_pl,Rd = A_a*f_yd + A_c*f_cd + A_s*f_sd")
print(f"N_pl,Rd = {Npl_Rd_a/1000:.1f} kN + {Npl_Rd_c/1000:.1f} kN + {Npl_Rd_s/1000:.1f} kN = {N_pl_Rd/1000:.2f} kN\n")

# Step 5: Calculate Effective Flexural Stiffness ((EI)_eff)
print("Step 5: Effective Flexural Stiffness (about weak z-z axis)")
# Secant modulus of elasticity of concrete
E_cm = 22000 * ((f_ck + 8) / 10)**0.3
# Moment of inertia of concrete (gross section) about z-z axis
I_cz = h_c * b_c**3 / 12
# Moment of inertia of reinforcement about z-z axis
d_rebar_pos = b_c/2 - cover - d_rebar/2
I_s = num_rebars * A_s_one_rebar * d_rebar_pos**2

# Effective flexural stiffness (EN 1994-1-1 eq. 6.30)
EI_eff = E_a * I_az + E_s * I_s + K_e * E_cm * I_cz
print(f"E_cm = {E_cm:.2f} N/mm^2")
print(f"I_az = {I_az:.2e} mm^4")
print(f"I_s = {I_s:.2e} mm^4")
print(f"I_cz = {I_cz:.2e} mm^4")
print(f"(EI)_eff = E_a*I_az + E_s*I_s + K_e*E_cm*I_cz = {EI_eff:.2e} N.mm^2\n")

# Step 6: Calculate Critical Buckling Load (N_cr) and Slenderness
print("Step 6: Critical Load and Slenderness")
N_cr = (math.pi**2 * EI_eff) / (L_cr**2)
print(f"N_cr = pi^2 * (EI)_eff / L_cr^2 = {N_cr/1000:.2f} kN")

# Characteristic plastic resistance (for slenderness calculation)
N_pl_Rk_a = A_a * f_y
N_pl_Rk_c = alpha_cc * A_c * f_ck
N_pl_Rk_s = A_s * f_sk
N_pl_Rk = N_pl_Rk_a + N_pl_Rk_c + N_pl_Rk_s
print(f"N_pl,Rk = {N_pl_Rk/1000:.2f} kN")

# Non-dimensional slenderness (EN 1994-1-1 eq. 6.29)
lambda_bar = math.sqrt(N_pl_Rk / N_cr)
print(f"Non-dimensional slenderness λ_bar = sqrt(N_pl,Rk / N_cr) = {lambda_bar:.3f}\n")

# Step 7: Calculate Reduction Factor (chi)
print("Step 7: Buckling Reduction Factor")
# For concrete encased I-section, weak axis (z-z), use buckling curve 'b' (EN 1994-1-1 Table 6.5)
# For curve 'b', imperfection factor alpha = 0.34 (EN 1993-1-1 Table 6.1)
alpha_imperfection = 0.34
# For composite columns, plateau length lambda_0 = 0.2
lambda_0 = 0.2

Phi = 0.5 * (1 + alpha_imperfection * (lambda_bar - lambda_0) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
print(f"Using Buckling Curve 'b' (alpha = {alpha_imperfection})")
print(f"Φ = 0.5 * [1 + α*(λ_bar - 0.2) + λ_bar^2] = {Phi:.3f}")
print(f"χ = 1 / (Φ + sqrt(Φ^2 - λ_bar^2)) = {chi:.3f}\n")

# Step 8: Calculate Design Buckling Resistance (N_b,Rd)
print("Step 8: Final Design Buckling Resistance")
N_b_Rd = chi * N_pl_Rd
print(f"The design buckling resistance is given by N_b,Rd = χ * N_pl,Rd")
print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd/1000:.2f} kN")
print(f"N_b,Rd = {N_b_Rd/1000:.2f} kN")

final_answer = N_b_Rd/1000
# "<<<" marks the start of the answer, and ">>>" marks the end.
print(f"\nFinal Answer:")
print(f"The buckling resistance of the column is {final_answer:.1f} kN.")
print(f"<<<{final_answer:.1f}>>>")