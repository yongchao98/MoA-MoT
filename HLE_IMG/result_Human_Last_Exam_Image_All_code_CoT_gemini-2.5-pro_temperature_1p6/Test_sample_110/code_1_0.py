import math

# Step 1: Define Material and Geometric Properties
# All units are in N and mm
print("Step 1: Define Material and Geometric Properties\n")

# Concrete C40/50
f_ck = 40.0  # Characteristic cylinder strength (MPa or N/mm^2)
gamma_c = 1.5  # Partial safety factor for concrete

# Structural Steel UC 254x254x132 S355
f_y = 345.0  # Yield strength (MPa or N/mm^2)
gamma_a = 1.0  # Partial safety factor for structural steel
A_a = 168.0 * 100  # Area from cm^2 to mm^2
I_ay = 22500.0 * 10**4  # Second moment of area y-y from cm^4 to mm^4
I_az = 7530.0 * 10**4   # Second moment of area z-z from cm^4 to mm^4
E_a = 210000.0  # Modulus of Elasticity for steel (N/mm^2)

# Reinforcement Steel S500
f_sk = 500.0  # Characteristic yield strength (N/mm^2)
gamma_s = 1.15  # Partial safety factor for rebar
d_rebar = 16.0  # Diameter of rebar (mm)
num_rebars = 4
E_s = 200000.0 # Modulus of Elasticity for rebar (N/mm^2)

# Column geometry
h_col = 400.0  # Column height (mm)
b_col = 400.0  # Column width (mm)
cover = 30.0  # Concrete cover (mm)
L_col = 4000.0  # Column effective length (mm)
L_cr = 1.0 * L_col # Buckling length (pinned-pinned assumption)

print(f"Concrete f_ck: {f_ck} N/mm^2")
print(f"Structural Steel f_y: {f_y} N/mm^2, Area: {A_a} mm^2")
print(f"Reinforcement f_sk: {f_sk} N/mm^2, Diameter: {d_rebar} mm")
print(f"Column Dimensions: {b_col} mm x {h_col} mm, Length: {L_cr} mm\n")

# Step 2: Calculate Design Strengths and Section Areas
print("Step 2: Calculate Design Strengths and Section Areas\n")

# Design strengths
f_cd = f_ck / gamma_c
f_yd = f_y / gamma_a
f_sd = f_sk / gamma_s
print(f"f_cd = {f_ck} / {gamma_c} = {f_cd:.2f} N/mm^2")
print(f"f_yd = {f_y} / {gamma_a} = {f_yd:.2f} N/mm^2")
print(f"f_sd = {f_sk} / {gamma_s} = {f_sd:.2f} N/mm^2\n")

# Areas
A_s_one = math.pi * (d_rebar / 2)**2
A_s = num_rebars * A_s_one
A_c_gross = h_col * b_col
A_c = A_c_gross - A_a - A_s
print(f"Area of Steel Section (A_a): {A_a:.2f} mm^2")
print(f"Area of Reinforcement (A_s): {num_rebars} * pi * ({d_rebar}/2)^2 = {A_s:.2f} mm^2")
print(f"Area of Concrete (A_c): {A_c_gross:.2f} - {A_a:.2f} - {A_s:.2f} = {A_c:.2f} mm^2\n")


# Step 3: Calculate Plastic Resistance of the Cross-Section (N_pl,Rd)
print("Step 3: Calculate Plastic Resistance (N_pl,Rd)\n")

# Contribution from each component
N_pl_Rd_a = A_a * f_yd
N_pl_Rd_c = 0.85 * A_c * f_cd  # 0.85 factor for concrete in encased sections
N_pl_Rd_s = A_s * f_sd
N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s

print("N_pl,Rd = A_a*f_yd + 0.85*A_c*f_cd + A_s*f_sd")
print(f"N_pl,Rd = {A_a:.0f}*{f_yd:.0f} + 0.85*{A_c:.0f}*{f_cd:.2f} + {A_s:.0f}*{f_sd:.2f}")
print(f"N_pl,Rd = {N_pl_Rd_a/1000:.1f} kN + {N_pl_Rd_c/1000:.1f} kN + {N_pl_Rd_s/1000:.1f} kN")
print(f"N_pl,Rd = {N_pl_Rd/1000:.1f} kN\n")

# Also calculate characteristic plastic resistance (N_pl,Rk) for slenderness calculation
N_pl_Rk = A_a * f_y + 0.85 * A_c * f_ck + A_s * f_sk
print(f"Characteristic Plastic Resistance (N_pl,Rk) = {N_pl_Rk/1000:.1f} kN\n")


# Step 4: Calculate Effective Flexural Stiffness ((EI)_eff)
print("Step 4: Calculate Effective Flexural Stiffness ((EI)_eff) for weaker z-z axis\n")

# Secant Modulus of Elasticity of Concrete
E_cm = 22 * ((f_ck + 8) / 10)**0.3
print(f"Concrete Modulus (E_cm) = 22 * (({f_ck} + 8) / 10)^0.3 = {E_cm:.0f} N/mm^2\n")

# Second moment of area of gross concrete section
I_c_gross_z = (h_col * b_col**3) / 12

# Second moment of area of reinforcement
# Distance of rebar centroid from column centroid
d_rebar_centroid = b_col / 2 - cover - d_rebar / 2
I_sz = num_rebars * A_s_one * d_rebar_centroid**2

# Effective Flexural Stiffness (z-z axis is weaker for UC section)
# (EI)_eff = E_a*I_a + K_e*E_cm*I_c + E_s*I_s
# K_e = 0.6 for short term loading
K_e = 0.6
EI_eff_z = E_a * I_az + K_e * E_cm * I_c_gross_z + E_s * I_sz
print("(EI)_eff,z = E_a*I_az + 0.6*E_cm*I_c,gross,z + E_s*I_sz")
print(f"(EI)_eff,z = {E_a:.0f}*{I_az:.0f} + {K_e}*{E_cm:.0f}*{I_c_gross_z:.0f} + {E_s:.0f}*{I_sz:.0f}")
print(f"(EI)_eff,z = {EI_eff_z:.2e} N.mm^2\n")


# Step 5: Calculate Critical Buckling Load (N_cr)
print("Step 5: Calculate Critical Buckling Load (N_cr)\n")

N_cr = (math.pi**2 * EI_eff_z) / (L_cr**2)
print(f"N_cr = pi^2 * (EI)_eff / L_cr^2")
print(f"N_cr = pi^2 * {EI_eff_z:.2e} / {L_cr:.0f}^2 = {N_cr:.0f} N")
print(f"N_cr = {N_cr/1000:.1f} kN\n")


# Step 6: Calculate Relative Slenderness (lambda_bar)
print("Step 6: Calculate Relative Slenderness (lambda_bar)\n")

lambda_bar = math.sqrt(N_pl_Rk / N_cr)
print("lambda_bar = sqrt(N_pl,Rk / N_cr)")
print(f"lambda_bar = sqrt({N_pl_Rk:.0f} / {N_cr:.0f}) = {lambda_bar:.4f}\n")


# Step 7: Calculate Buckling Reduction Factor (chi)
print("Step 7: Calculate Buckling Reduction Factor (chi)\n")

# For encased I-section buckling about weak z-z axis, use buckling curve 'b'
alpha = 0.34 # Imperfection factor for buckling curve 'b'
print(f"Using buckling curve 'b', imperfection factor alpha = {alpha}\n")

phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
print("Phi = 0.5 * [1 + alpha * (lambda_bar - 0.2) + lambda_bar^2]")
print(f"Phi = 0.5 * [1 + {alpha} * ({lambda_bar:.4f} - 0.2) + {lambda_bar:.4f}^2] = {phi:.4f}")
print("chi = 1 / (Phi + sqrt(Phi^2 - lambda_bar^2))")
print(f"chi = 1 / ({phi:.4f} + sqrt({phi:.4f}^2 - {lambda_bar:.4f}^2)) = {chi:.4f}\n")


# Step 8: Calculate Design Buckling Resistance (N_b,Rd)
print("Step 8: Calculate Final Design Buckling Resistance (N_b,Rd)\n")

N_b_Rd = chi * N_pl_Rd
print("N_b,Rd = chi * N_pl,Rd")
print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd/1000:.1f} kN")
print(f"N_b,Rd = {N_b_Rd/1000:.1f} kN")


print("\n--- FINAL ANSWER ---")
print(f"The buckling resistance of the column is {N_b_Rd / 1000:.1f} kN.")
# Final answer in specified format
final_answer = N_b_Rd / 1000
print(f'<<<{final_answer:.1f}>>>')
