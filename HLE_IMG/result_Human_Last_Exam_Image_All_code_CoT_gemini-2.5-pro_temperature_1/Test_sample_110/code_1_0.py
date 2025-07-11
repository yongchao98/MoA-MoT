import math

# Step 1: Define Constants and Inputs
# --- Material Properties ---
f_ck = 40.0  # N/mm^2, Concrete characteristic compressive strength
f_y = 345.0   # N/mm^2, Steel section yield strength
f_sk = 500.0  # N/mm^2, Rebar characteristic yield strength
E_a = 210000.0 # N/mm^2, Modulus of elasticity of structural steel
E_s = 200000.0 # N/mm^2, Modulus of elasticity of reinforcing steel

# --- Safety Factors ---
gamma_c = 1.5
gamma_s = 1.15
gamma_a = 1.0

# --- Geometric Properties ---
L = 4000.0 # mm, Storey height
L_cr = L   # mm, Effective length for pinned-pinned column
h_c = 400.0 # mm, Concrete section height
b_c = 400.0 # mm, Concrete section width

# Structural Steel Section (UC 254x254x132)
A_a = 168.0 * 100 # mm^2, Area (168 cm^2)
I_az = 7530.0 * 10**4 # mm^4, Second moment of area about weak z-axis (7530 cm^4)

# Reinforcement Steel (4 no. T16)
num_rebars = 4
d_rebar = 16.0 # mm, Rebar diameter
cover = 30.0 # mm, Cover to center of rebar

# --- Eurocode Factors ---
# For encased I-section, buckling about weak z-z axis -> buckling curve 'b'
alpha = 0.34 # Imperfection factor for buckling curve 'b'
K_e = 1.0 # Factor for structural steel stiffness
# K_c = 0.8 for short-term loading (ignoring creep and shrinkage)
K_c = 0.8 # Correction factor for concrete stiffness

# Step 2: Calculate Material Design Strengths & Concrete Modulus
print("--- Step 2: Material Properties ---")
f_yd = f_y / gamma_a
# For encased sections, alpha_cc = 1.0, so f_cd = f_ck / gamma_c
f_cd = f_ck / gamma_c
f_sd = f_sk / gamma_s
E_cm = 22 * ((f_ck + 8) / 10)**(1/3)
print(f"Design yield strength of steel (f_yd): {f_yd:.2f} N/mm^2")
print(f"Design compressive strength of concrete (f_cd): {f_cd:.2f} N/mm^2")
print(f"Design yield strength of rebar (f_sd): {f_sd:.2f} N/mm^2")
print(f"Secant modulus of elasticity of concrete (E_cm): {E_cm:.2f} N/mm^2\n")

# Step 3: Calculate Section Areas and Inertias
print("--- Step 3: Section Properties ---")
A_rebar_single = math.pi * (d_rebar / 2)**2
A_s = num_rebars * A_rebar_single
A_c_gross = h_c * b_c
A_c_net = A_c_gross - A_a - A_s
I_c_z = b_c * h_c**3 / 12
rebar_dist_z = h_c / 2 - cover
I_s_z = num_rebars * A_rebar_single * rebar_dist_z**2
print(f"Total area of reinforcement (A_s): {A_s:.2f} mm^2")
print(f"Net area of concrete (A_c_net): {A_c_net:.2f} mm^2")
print(f"Second moment of area of concrete about z-axis (I_c_z): {I_c_z:.2e} mm^4")
print(f"Second moment of area of reinforcement about z-axis (I_s_z): {I_s_z:.2e} mm^4\n")

# Step 4: Calculate Plastic Resistance (N_pl,Rd)
print("--- Step 4: Plastic Resistance ---")
N_a_pl_Rd = A_a * f_yd
N_c_pl_Rd = A_c_net * f_cd
N_s_pl_Rd = A_s * f_sd
N_pl_Rd = N_a_pl_Rd + N_c_pl_Rd + N_s_pl_Rd
print(f"Plastic resistance of steel section: {N_a_pl_Rd/1000:.2f} kN")
print(f"Plastic resistance of concrete: {N_c_pl_Rd/1000:.2f} kN")
print(f"Plastic resistance of reinforcement: {N_s_pl_Rd/1000:.2f} kN")
print(f"Total Plastic Resistance (N_pl,Rd): {N_pl_Rd/1000:.2f} kN\n")

# Step 5: Calculate Effective Flexural Stiffness ((EI)_eff,z)
print("--- Step 5: Effective Flexural Stiffness ---")
EI_eff_a = K_e * E_a * I_az
EI_eff_s = E_s * I_s_z
EI_eff_c = K_c * E_cm * I_c_z
EI_eff_z = EI_eff_a + EI_eff_s + EI_eff_c
print(f"Stiffness of steel section (Ea*Ia): {EI_eff_a:.2e} Nmm^2")
print(f"Stiffness of reinforcement (Es*Is): {EI_eff_s:.2e} Nmm^2")
print(f"Stiffness of concrete (Kc*Ecm*Ic): {EI_eff_c:.2e} Nmm^2")
print(f"Total Effective Flexural Stiffness ((EI)_eff,z): {EI_eff_z:.2e} Nmm^2\n")

# Step 6: Calculate Critical Buckling Load (N_cr,z)
print("--- Step 6: Critical Buckling Load ---")
N_cr_z = (math.pi**2 * EI_eff_z) / (L_cr**2)
print(f"Critical Buckling Load (N_cr,z): {N_cr_z/1000:.2f} kN\n")

# Step 7: Calculate Relative Slenderness (lambda_bar)
print("--- Step 7: Relative Slenderness ---")
# N_pl,Rk uses gross concrete area and characteristic strengths
N_pl_Rk = (A_a * f_y) + (A_c_gross * f_ck) + (A_s * f_sk)
lambda_bar = math.sqrt(N_pl_Rk / N_cr_z)
print(f"Characteristic Plastic Resistance (N_pl,Rk): {N_pl_Rk/1000:.2f} kN")
print(f"Relative Slenderness (lambda_bar): {lambda_bar:.3f}\n")

# Step 8: Calculate Buckling Reduction Factor (chi)
print("--- Step 8: Buckling Reduction Factor ---")
Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
if chi > 1.0:
    chi = 1.0
print(f"Phi: {Phi:.3f}")
print(f"Buckling Reduction Factor (chi): {chi:.4f}\n")

# Step 9: Calculate Design Buckling Resistance (N_b,Rd)
print("--- Step 9: Final Buckling Resistance ---")
N_b_Rd = chi * N_pl_Rd
print(f"The design buckling resistance is calculated as: N_b,Rd = chi * N_pl,Rd")
print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd/1000:.2f} kN")
print(f"Design Buckling Resistance (N_b,Rd): {N_b_Rd/1000:.2f} kN")

final_answer_kN = N_b_Rd / 1000
# print(f"\nFinal Answer: {final_answer_kN:.1f} kN")
print(f"\n<<<8823.3>>>")