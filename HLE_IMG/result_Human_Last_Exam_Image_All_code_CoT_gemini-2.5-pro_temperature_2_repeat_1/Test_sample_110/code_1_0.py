import math

# Step 1: Define Material Properties
# ---
print("Step 1: Material Properties")
# Concrete C40/50
f_ck = 40.0  # N/mm^2, characteristic cylinder strength
gamma_c = 1.5  # Partial factor for concrete
alpha_cc = 0.85 # Coefficient for long-term effects
f_cd = alpha_cc * f_ck / gamma_c
E_cm = 22 * ((f_ck + 8) / 10)**0.3 # N/mm^2, secant modulus of elasticity of concrete

# Structural Steel UC 254x254x132 S355
f_y = 345.0  # N/mm^2, yield strength
gamma_a = 1.0  # Partial factor for structural steel
f_yd = f_y / gamma_a
E_a = 210000.0  # N/mm^2, modulus of elasticity

# Reinforcement Steel S500
f_sk = 500.0  # N/mm^2, characteristic yield strength
gamma_s = 1.15 # Partial factor for reinforcing steel
f_sd = f_sk / gamma_s
E_s = 200000.0  # N/mm^2, modulus of elasticity

print(f"Design strength of concrete f_cd = {alpha_cc:.2f} * {f_ck:.1f} / {gamma_c:.2f} = {f_cd:.2f} N/mm^2")
print(f"Modulus of elasticity of concrete E_cm = 22 * (({f_ck:.1f} + 8) / 10)^0.3 = {E_cm:.0f} N/mm^2")
print(f"Design strength of structural steel f_yd = {f_y:.1f} / {gamma_a:.1f} = {f_yd:.1f} N/mm^2")
print(f"Design strength of reinforcement f_sd = {f_sk:.1f} / {gamma_s:.2f} = {f_sd:.2f} N/mm^2\n")

# Step 2: Define Sectional Properties (weaker z-z axis)
# ---
print("Step 2: Sectional Properties (for weaker z-z axis)")
# Concrete Section
h_c = 400.0  # mm
b_c = 400.0  # mm
A_c = h_c * b_c  # mm^2, gross area of concrete
I_cz = b_c * h_c**3 / 12  # mm^4, second moment of area of concrete about z-axis

# Structural Steel Profile
A_a = 168.0 * 100 # mm^2, Area from cm^2
I_az = 7530.0 * 10000 # mm^4, Second moment of area from cm^4

# Reinforcement Steel (4 x T16)
d_bar = 16.0 # mm, bar diameter
cover = 30.0 # mm
A_s_one_bar = math.pi * (d_bar / 2)**2
A_s = 4 * A_s_one_bar # mm^2, total area of reinforcement
# Distance from section centroid (z-axis) to rebar centroid
y_s = h_c/2 - cover - d_bar/2
# Second moment of area of reinforcement about z-axis (Parallel Axis Theorem)
I_sz = A_s * y_s**2

print(f"Gross area of concrete Ac = {h_c:.1f} mm * {b_c:.1f} mm = {A_c:.0f} mm^2")
print(f"Total area of reinforcement As = 4 * pi * ({d_bar:.1f}/2)^2 = {A_s:.1f} mm^2")
print(f"Second moment of area of concrete I_cz = {b_c:.1f}*({h_c:.1f}^3)/12 = {I_cz:.0f} mm^4")
print(f"Second moment of area of steel I_az = {I_az:.0f} mm^4")
print(f"Distance from centroid to rebar y_s = {h_c/2:.1f} - {cover:.1f} - {d_bar/2:.1f} = {y_s:.1f} mm")
print(f"Second moment of area of rebar I_sz = {A_s:.1f} mm^2 * ({y_s:.1f} mm)^2 = {I_sz:.0f} mm^4\n")

# Step 3: Calculate Effective Flexural Stiffness (EI)_eff,z
# ---
print("Step 3: Calculate Effective Flexural Stiffness (EI)_eff,z")
# According to EN 1994-1-1 6.7.3.3(3), for encased sections, a reduction factor k_c = 0.6 is applied to I_c.
k_c = 0.6
EI_eff_z = E_a * I_az + k_c * E_cm * I_cz + E_s * I_sz

print(f"(EI)_eff,z = E_a*I_az + k_c*E_cm*I_cz + E_s*I_sz")
print(f"(EI)_eff,z = ({E_a/1e3:.0f} * {I_az/1e6:.1f}) + ({k_c:.1f} * {E_cm/1e3:.1f} * {I_cz/1e6:.1f}) + ({E_s/1e3:.0f} * {I_sz/1e6:.1f}) [kNm^2]")
print(f"(EI)_eff,z = {(E_a*I_az)/1e9:.1f} + {(k_c*E_cm*I_cz)/1e9:.1f} + {(E_s*I_sz)/1e9:.1f} [kNm^2]")
print(f"(EI)_eff,z = {EI_eff_z / 1e9:.1f} kNm^2\n")

# Step 4: Calculate Critical Buckling Load N_cr,z
# ---
print("Step 4: Calculate Critical Buckling Load N_cr,z")
L = 4.0 # m, storey height
k_buckling = 1.0 # Effective length factor for pinned-pinned column
L_cr = k_buckling * L * 1000 # mm
N_cr_z = (math.pi**2 * EI_eff_z) / L_cr**2

print(f"Effective buckling length L_cr = {k_buckling:.1f} * {L*1000:.0f} mm = {L_cr:.0f} mm")
print(f"N_cr,z = (pi^2 * (EI)_eff,z) / L_cr^2 = (pi^2 * {EI_eff_z/1e9:.1f} kNm^2) / ({L_cr/1000:.1f} m)^2 = {N_cr_z / 1e3:.1f} kN\n")

# Step 5: Calculate Plastic Resistances N_pl,Rk and N_pl,Rd
# ---
print("Step 5: Calculate Plastic Resistances")
# Per EN1994-1-1 6.7.3.2(1) note, for encased sections a factor of 0.85 is used for concrete contribution.
# Design plastic resistance N_pl,Rd
N_pl_Rd = A_a * f_yd + 0.85 * A_c * f_cd + A_s * f_sd
# Characteristic plastic resistance N_pl,Rk
N_pl_Rk = A_a * f_y + 0.85 * A_c * f_ck + A_s * f_sk

print(f"N_pl,Rd = Aa*f_yd + 0.85*Ac*f_cd + As*f_sd")
print(f"N_pl,Rd = ({A_a:.0f}*{f_yd:.1f}) + (0.85*{A_c:.0f}*{f_cd:.2f}) + ({A_s:.1f}*{f_sd:.2f})")
print(f"N_pl,Rd = {(A_a * f_yd)/1e3:.0f} + {(0.85 * A_c * f_cd)/1e3:.0f} + {(A_s * f_sd)/1e3:.0f} = {N_pl_Rd / 1e3:.1f} kN")
print(f"N_pl,Rk = Aa*f_y + 0.85*Ac*f_ck + As*f_sk = {N_pl_Rk / 1e3:.1f} kN\n")


# Step 6: Determine Buckling Reduction Factor chi
# ---
print("Step 6: Determine Buckling Reduction Factor chi")
# Relative slenderness lambda_bar
lambda_bar = math.sqrt(N_pl_Rk / N_cr_z)

# For buckling about z-z axis of an encased I-section, use Buckling Curve 'b' (EN1994-1-1 Table 6.5)
# Imperfection factor for curve 'b' is alpha = 0.34 (EN1993-1-1 Table 6.2)
alpha = 0.34

# Intermediate factor phi
phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)

# Buckling reduction factor chi
chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))

print(f"Relative Slenderness lambda_bar = sqrt(N_pl,Rk / N_cr,z) = sqrt({N_pl_Rk/1e3:.1f} / {N_cr_z/1e3:.1f}) = {lambda_bar:.4f}")
print(f"For buckling curve 'b', alpha = {alpha}")
print(f"phi = 0.5 * [1 + alpha*(lambda_bar - 0.2) + lambda_bar^2] = 0.5 * [1 + {alpha:.2f}*({lambda_bar:.4f} - 0.2) + {lambda_bar**2:.4f}] = {phi:.4f}")
print(f"Reduction factor chi = 1 / (phi + sqrt(phi^2 - lambda_bar^2)) = 1 / ({phi:.4f} + sqrt({phi**2:.4f} - {lambda_bar**2:.4f})) = {chi:.4f}\n")


# Step 7: Calculate Design Buckling Resistance N_b,Rd
# ---
print("Step 7: Final Buckling Resistance N_b,Rd")
N_b_Rd = chi * N_pl_Rd
print(f"N_b,Rd = chi * N_pl,Rd = {chi:.4f} * {N_pl_Rd / 1e3:.1f} kN = {N_b_Rd / 1e3:.1f} kN")

final_answer_kN = N_b_Rd / 1000
# print(f"\nFinal Answer: {final_answer_kN:.1f} kN")

<<<7998.7>>>