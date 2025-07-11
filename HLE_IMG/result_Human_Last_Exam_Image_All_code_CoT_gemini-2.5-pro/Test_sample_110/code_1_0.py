import math

# Step 1: Define Material and Section Properties from the problem statement
# --- Geometry
h_col = 400.0  # mm, height of concrete section
b_col = 400.0  # mm, width of concrete section
L_col = 4000.0 # mm, column height
cover = 30.0   # mm, concrete cover
d_rebar = 16.0 # mm, diameter of rebar
n_rebar = 4    # number of rebars

# --- Material Properties (Characteristic Strengths)
f_ck = 40.0  # MPa, concrete characteristic compressive strength
f_y = 345.0  # MPa, structural steel yield strength
f_sk = 500.0 # MPa, reinforcing steel yield strength

# --- Material Moduli of Elasticity
E_a = 210000.0 # MPa, Modulus of elasticity for structural steel
E_s = 210000.0 # MPa, Modulus of elasticity for reinforcing steel
E_cm = 35000.0 # MPa, Secant modulus of elasticity of concrete (C40/50)

# --- Structural Steel Section Properties (UC 254 x 254 x 132)
A_a = 168.0 * 100 # mm^2, Area (from 168 cm^2)
I_ay = 22500.0 * 10**4 # mm^4, Second moment of area about y-axis (from 22500 cm^4)
I_az = 7530.0 * 10**4  # mm^4, Second moment of area about z-axis (from 7530 cm^4)

# --- Partial Safety Factors (Eurocode)
gamma_a = 1.0
gamma_c = 1.5
gamma_s = 1.15

# --- Other Factors
K_e = 0.6  # Factor for concrete stiffness (ignoring long-term effects)
k_buckling = 1.0 # Effective length factor (pinned-pinned ends assumed)

# Step 2: Calculate Design Strengths and Areas
f_yd = f_y / gamma_a
f_cd = f_ck / gamma_c
f_sd = f_sk / gamma_s

A_rebar_single = math.pi * (d_rebar / 2)**2
A_s = n_rebar * A_rebar_single
A_gross = h_col * b_col
A_c = A_gross - A_a - A_s

# Step 3: Calculate Plastic Resistance (N_pl,Rd)
N_pl_Rd_a = A_a * f_yd
N_pl_Rd_c = 0.85 * A_c * f_cd
N_pl_Rd_s = A_s * f_sd
N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s

# Step 4: Calculate Effective Flexural Stiffness ((EI)_eff)
# Buckling occurs about the weaker z-z axis (since I_az < I_ay)
I_c = (b_col * h_col**3) / 12
# Distance of rebar center from the section's centroidal z-axis
d_rebar_z = h_col / 2 - cover - d_rebar / 2
I_s = n_rebar * A_rebar_single * d_rebar_z**2 # Using parallel axis theorem, ignoring rebar's own I

EI_eff_a = E_a * I_az
EI_eff_s = E_s * I_s
EI_eff_c = K_e * E_cm * I_c
EI_eff = EI_eff_a + EI_eff_s + EI_eff_c

# Step 5: Calculate Critical Buckling Load (N_cr)
L_cr = k_buckling * L_col
N_cr = (math.pi**2 * EI_eff) / L_cr**2

# Step 6: Determine Relative Slenderness (λ_bar)
# First, calculate characteristic plastic resistance N_pl,Rk
N_pl_Rk = A_a * f_y + A_c * f_ck + A_s * f_sk
lambda_bar = math.sqrt(N_pl_Rk / N_cr)

# Step 7: Calculate Buckling Reduction Factor (χ)
# For encased I-section, weak axis (z-z), use buckling curve 'c' -> alpha = 0.49
alpha = 0.49
Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
# Ensure chi is not greater than 1.0
chi = min(chi, 1.0)

# Step 8: Calculate Final Buckling Resistance (N_b,Rd)
N_b_Rd = chi * N_pl_Rd

# --- Output the results ---
print("--- Calculation of Buckling Resistance ---\n")
print(f"1. Plastic Resistance (N_pl,Rd)")
print(f"   N_pl,Rd = {N_pl_Rd / 1000:.1f} kN\n")

print(f"2. Effective Flexural Stiffness ((EI)_eff)")
print(f"   (EI)_eff = {EI_eff:.3e} N.mm^2\n")

print(f"3. Critical Buckling Load (N_cr)")
print(f"   N_cr = {N_cr / 1000:.1f} kN\n")

print(f"4. Relative Slenderness (λ_bar)")
print(f"   λ_bar = {lambda_bar:.3f}\n")

print(f"5. Buckling Reduction Factor (χ)")
print(f"   Φ = {Phi:.3f}")
print(f"   χ = {chi:.3f}\n")

print("--- Final Buckling Resistance ---\n")
print("The buckling resistance is calculated as:")
print("N_b,Rd = χ * N_pl,Rd")
print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd/1000:.1f} kN")
print(f"N_b,Rd = {N_b_Rd/1000:.1f} kN")
print("\nThe buckling resistance of the column is:")
print(f"{N_b_Rd/1000:.1f} kN")
print(f"\n<<<7667.6>>>")