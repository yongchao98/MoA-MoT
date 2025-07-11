import math

# Step 1: Define Material and Sectional Properties
# All units are in N and mm

# Steel Section: UC 254 x 254 x 132 S355
fy = 345.0      # N/mm^2 (Yield strength)
A_a = 168.0e2   # mm^2 (Area, 168 cm^2)
I_ay = 22500.0e4# mm^4 (Moment of inertia about y-axis, 22500 cm^4)
I_az = 7530.0e4 # mm^4 (Moment of inertia about z-axis, 7530 cm^4)
E_a = 210000.0  # N/mm^2 (Young's Modulus for structural steel)

# Reinforcing Steel: 4 no. T16 S500
fsk = 500.0     # N/mm^2 (Characteristic yield strength of reinforcement)
n_s = 4         # Number of rebars
d_s = 16.0      # mm (Diameter of rebar)
A_s1 = math.pi * (d_s / 2)**2 # mm^2 (Area of one rebar)
A_s = n_s * A_s1# mm^2 (Total area of reinforcement)
E_s = 200000.0  # N/mm^2 (Young's Modulus for reinforcing steel)

# Concrete: C40/50
fck = 40.0      # N/mm^2 (Characteristic cylinder strength)
E_cm = 22 * ((fck + 8) / 10)**0.3 # N/mm^2 (Secant modulus of elasticity)

# Section Geometry and other parameters
h_c = 400.0     # mm (Overall height of section)
b_c = 400.0     # mm (Overall width of section)
cover = 30.0    # mm (Concrete cover)
L = 4000.0      # mm (Storey height / Column length)
L_cr = L * 1.0  # mm (Effective buckling length, assuming pinned-pinned k=1.0)
K_e = 0.6       # Correction factor for concrete modulus
alpha_cc = 0.85 # Factor for compressive strength of concrete in composite sections

# Partial safety factors
gamma_a = 1.0   # Structural steel
gamma_s = 1.15  # Reinforcing steel
gamma_c = 1.5   # Concrete

# Step 2: Calculate Plastic Resistance of Cross-section (N_pl,Rd)
A_c = h_c * b_c - A_a - A_s # mm^2 (Net area of concrete)
f_yd = fy / gamma_a
f_sd = fsk / gamma_s
f_cd = fck / gamma_c

N_pl_Rd = A_a * f_yd + alpha_cc * A_c * f_cd + A_s * f_sd
print(f"--- Calculating Plastic Resistance (N_pl,Rd) ---")
print(f"Area of Steel section (A_a): {A_a:.2f} mm^2")
print(f"Area of Rebar (A_s): {A_s:.2f} mm^2")
print(f"Area of Concrete (A_c): {A_c:.2f} mm^2")
print(f"Design yield strength of steel (f_yd): {f_yd:.2f} N/mm^2")
print(f"Design strength of concrete (f_cd): {f_cd:.2f} N/mm^2")
print(f"Design strength of rebar (f_sd): {f_sd:.2f} N/mm^2")
print(f"Plastic Resistance N_pl,Rd = (A_a * f_yd) + (alpha_cc * A_c * f_cd) + (A_s * f_sd)")
print(f"N_pl,Rd = ({A_a:.0f} * {f_yd:.0f}) + ({alpha_cc} * {A_c:.0f} * {f_cd:.2f}) + ({A_s:.0f} * {f_sd:.2f})")
print(f"N_pl,Rd = {A_a * f_yd:.0f} + {alpha_cc * A_c * f_cd:.0f} + {A_s * f_sd:.0f} = {N_pl_Rd:.0f} N")
print(f"Plastic Resistance N_pl,Rd = {N_pl_Rd / 1000:.1f} kN\n")

# Step 3: Calculate Effective Flexural Stiffness (EI)_eff
# Buckling will occur about the weaker axis (z-z)
I_c = b_c * h_c**3 / 12  # Gross moment of inertia of concrete section
# Distance from centroid to rebar center
d_rebar_z = h_c/2 - cover - d_s/2
I_s = n_s * A_s1 * d_rebar_z**2 # Moment of inertia of reinforcement

EI_eff_z = E_a * I_az + K_e * E_cm * I_c + E_s * I_s
print(f"--- Calculating Effective Flexural Stiffness ((EI)_eff,z) ---")
print(f"Moment of inertia of steel section (I_az): {I_az:.2e} mm^4")
print(f"Moment of inertia of concrete (I_c): {I_c:.2e} mm^4")
print(f"Moment of inertia of rebar (I_s): {I_s:.2e} mm^4")
print(f"(EI)_eff,z = E_a*I_az + K_e*E_cm*I_c + E_s*I_s = {EI_eff_z:.2e} N.mm^2\n")

# Step 4: Calculate Elastic Critical Buckling Load (N_cr)
N_cr_z = (math.pi**2 * EI_eff_z) / L_cr**2
print(f"--- Calculating Elastic Critical Buckling Load (N_cr) ---")
print(f"Effective length (L_cr): {L_cr:.0f} mm")
print(f"N_cr = (pi^2 * (EI)_eff) / L_cr^2 = (pi^2 * {EI_eff_z:.2e}) / {L_cr:.0f}^2 = {N_cr_z / 1000:.1f} kN\n")

# Step 5: Calculate Non-dimensional Slenderness (lambda_bar)
# N_pl_Rk is the characteristic plastic resistance
N_pl_Rk = A_a * fy + alpha_cc * A_c * fck + A_s * fsk
lambda_bar = math.sqrt(N_pl_Rk / N_cr_z)
print(f"--- Calculating Non-dimensional Slenderness (lambda_bar) ---")
print(f"Characteristic Plastic Resistance (N_pl,Rk): {N_pl_Rk / 1000:.1f} kN")
print(f"Lambda_bar = sqrt(N_pl_Rk / N_cr) = sqrt({N_pl_Rk:.0f} / {N_cr_z:.0f}) = {lambda_bar:.4f}\n")

# Step 6: Determine Buckling Reduction Factor (chi)
# For encased I-section buckling about weak axis z-z, use buckling curve 'b'
alpha_imp = 0.34 # Imperfection factor for curve 'b'
Phi = 0.5 * (1 + alpha_imp * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
print(f"--- Calculating Buckling Reduction Factor (chi) ---")
print(f"Using buckling curve 'b', imperfection factor (alpha) = {alpha_imp}")
print(f"Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar^2) = {Phi:.4f}")
print(f"chi = 1 / (Phi + sqrt(Phi^2 - lambda_bar^2)) = {chi:.4f}\n")

# Step 7: Calculate Buckling Resistance (N_b,Rd)
N_b_Rd = chi * N_pl_Rd
print(f"--- Final Calculation: Buckling Resistance (N_b,Rd) ---")
print("The buckling resistance is given by the equation: N_b,Rd = chi * N_pl,Rd")
print(f"chi = {chi:.4f}")
print(f"N_pl,Rd = {N_pl_Rd / 1000:.1f} kN")
print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd / 1000:.1f} kN")
print(f"\nThe buckling resistance of the column is {N_b_Rd / 1000:.1f} kN.")

final_answer = N_b_Rd / 1000
#<<<8177.6>>>