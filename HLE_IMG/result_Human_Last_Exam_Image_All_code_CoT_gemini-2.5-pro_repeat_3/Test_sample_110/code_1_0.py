import math

# Step 1: Define Material and Sectional Properties
# Given values from the problem description
# Column properties
storey_height = 4000.0  # mm
cover = 30.0  # mm
col_h = 400.0  # mm
col_b = 400.0  # mm

# Structural Steel (UC 254x254x132, S355)
f_y = 345.0  # N/mm^2
A_a = 168.0 * 100  # mm^2 (from 168 cm^2)
I_az = 7530.0 * 1e4  # mm^4 (from 7530 cm^4), weaker axis
E_a = 210000.0  # N/mm^2 (Modulus of Elasticity)

# Concrete (C40/50)
f_ck = 40.0  # N/mm^2

# Reinforcing Steel (4 no. T16, S500)
num_rebars = 4
rebar_dia = 16.0  # mm
f_sk = 500.0  # N/mm^2 (Characteristic yield strength)
E_s = 200000.0  # N/mm^2 (Modulus of Elasticity)

# Partial safety factors (from Eurocode)
gamma_a = 1.0
gamma_c = 1.5
gamma_s = 1.15

print("--- Step 1: Input Properties ---")
print(f"Storey Height (L): {storey_height} mm")
print(f"Steel Yield Strength (f_y): {f_y} N/mm^2")
print(f"Steel Area (A_a): {A_a} mm^2")
print(f"Steel Moment of Inertia (I_az): {I_az:.2e} mm^4")
print(f"Concrete Strength (f_ck): {f_ck} N/mm^2")
print(f"Rebar Strength (f_sk): {f_sk} N/mm^2")
print(f"Number of rebars: {num_rebars}, Diameter: {rebar_dia} mm")
print("-" * 30 + "\n")


# Step 2: Calculate Derived Material and Sectional Properties
# Design strengths
f_yd = f_y / gamma_a
f_cd = f_ck / gamma_c
f_sd = f_sk / gamma_s

# Concrete Modulus of Elasticity (E_cm)
E_cm = 22 * ((f_ck + 8) / 10)**0.3

# Concrete section properties
A_c_gross = col_h * col_b
A_c = A_c_gross - A_a
I_c_gross = (col_b * col_h**3) / 12

# Reinforcing steel properties
A_s1 = math.pi * (rebar_dia / 2)**2
A_s = num_rebars * A_s1
# Distance from section centroid to rebar centroid
d_s = col_h / 2 - cover - rebar_dia / 2
# Second moment of area of rebars about z-axis (using parallel axis theorem)
I_s = num_rebars * (A_s1 * d_s**2) # Self-inertia of bars is negligible

print("--- Step 2: Derived Properties ---")
print(f"Design Steel Strength (f_yd): {f_yd:.2f} N/mm^2")
print(f"Design Concrete Strength (f_cd): {f_cd:.2f} N/mm^2")
print(f"Design Rebar Strength (f_sd): {f_sd:.2f} N/mm^2")
print(f"Concrete Modulus of Elasticity (E_cm): {E_cm:.2f} N/mm^2")
print(f"Net Concrete Area (A_c): {A_c:.2f} mm^2")
print(f"Total Rebar Area (A_s): {A_s:.2f} mm^2")
print(f"Gross Concrete Moment of Inertia (I_c_gross): {I_c_gross:.2e} mm^4")
print(f"Rebar Moment of Inertia (I_s): {I_s:.2e} mm^4")
print("-" * 30 + "\n")


# Step 3: Calculate Plastic Resistance (N_pl,Rd)
N_pl_Rd = (A_a * f_yd) + (0.85 * A_c * f_cd) + (A_s * f_sd)
N_pl_Rk = (A_a * f_y) + (0.85 * A_c * f_ck) + (A_s * f_sk) # Characteristic resistance

print("--- Step 3: Plastic Resistance Calculation ---")
print(f"N_pl,Rd = A_a * f_yd + 0.85 * A_c * f_cd + A_s * f_sd")
print(f"N_pl,Rd = ({A_a} * {f_yd:.2f}) + (0.85 * {A_c:.2f} * {f_cd:.2f}) + ({A_s:.2f} * {f_sd:.2f})")
print(f"N_pl,Rd = {A_a * f_yd:.0f} + {0.85 * A_c * f_cd:.0f} + {A_s * f_sd:.0f} = {N_pl_Rd:.0f} N")
print(f"Design Plastic Resistance (N_pl,Rd): {N_pl_Rd / 1000:.2f} kN")
print(f"Characteristic Plastic Resistance (N_pl,Rk): {N_pl_Rk / 1000:.2f} kN")
print("-" * 30 + "\n")


# Step 4: Calculate Effective Flexural Stiffness ((EI)_eff)
# Using simplified method from EN1994-1-1, with K_e = 0.6 for concrete
K_e = 0.6
EI_eff = (E_a * I_az) + (K_e * E_cm * I_c_gross) + (E_s * I_s)

print("--- Step 4: Effective Flexural Stiffness Calculation ---")
print(f"(EI)_eff = E_a*I_az + K_e*E_cm*I_c + E_s*I_s")
print(f"(EI)_eff = ({E_a:.0f} * {I_az:.2e}) + ({K_e} * {E_cm:.2f} * {I_c_gross:.2e}) + ({E_s:.0f} * {I_s:.2e})")
print(f"(EI)_eff = {(E_a * I_az):.2e} + {(K_e * E_cm * I_c_gross):.2e} + {(E_s * I_s):.2e} = {EI_eff:.2e} Nmm^2")
print("-" * 30 + "\n")


# Step 5: Calculate Critical Elastic Buckling Load (N_cr)
# Assuming pinned-pinned ends, effective length L_cr = storey_height
L_cr = storey_height
N_cr = (math.pi**2 * EI_eff) / (L_cr**2)

print("--- Step 5: Critical Buckling Load Calculation ---")
print(f"N_cr = pi^2 * (EI)_eff / L_cr^2")
print(f"N_cr = ({math.pi**2:.2f} * {EI_eff:.2e}) / {L_cr}^2 = {N_cr:.0f} N")
print(f"Critical Buckling Load (N_cr): {N_cr / 1000:.2f} kN")
print("-" * 30 + "\n")


# Step 6: Determine Relative Slenderness (lambda_bar)
lambda_bar = math.sqrt(N_pl_Rk / N_cr)

print("--- Step 6: Relative Slenderness Calculation ---")
print(f"lambda_bar = sqrt(N_pl,Rk / N_cr)")
print(f"lambda_bar = sqrt({N_pl_Rk:.0f} / {N_cr:.0f}) = {lambda_bar:.3f}")
print("-" * 30 + "\n")


# Step 7: Find Buckling Reduction Factor (chi)
# For encased I-section, minor axis buckling -> use buckling curve 'c'
alpha = 0.49  # Imperfection factor for curve 'c'
Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
if chi > 1.0:
    chi = 1.0

print("--- Step 7: Buckling Reduction Factor Calculation ---")
print(f"Using Buckling Curve 'c' (alpha = {alpha})")
print(f"Phi = 0.5 * [1 + alpha * (lambda_bar - 0.2) + lambda_bar^2]")
print(f"Phi = 0.5 * [1 + {alpha} * ({lambda_bar:.3f} - 0.2) + {lambda_bar:.3f}^2] = {Phi:.3f}")
print(f"chi = 1 / (Phi + sqrt(Phi^2 - lambda_bar^2))")
print(f"chi = 1 / ({Phi:.3f} + sqrt({Phi:.3f}^2 - {lambda_bar:.3f}^2)) = {chi:.3f}")
print("-" * 30 + "\n")


# Step 8: Calculate Buckling Resistance (N_b,Rd)
N_b_Rd = chi * N_pl_Rd

print("--- Step 8: Final Buckling Resistance Calculation ---")
print(f"N_b,Rd = chi * N_pl,Rd")
print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd / 1000:.2f} kN")
final_answer_kN = N_b_Rd / 1000
print(f"\nThe final buckling resistance of the column is: {final_answer_kN:.2f} kN")
print(f"<<<{final_answer_kN:.2f}>>>")
