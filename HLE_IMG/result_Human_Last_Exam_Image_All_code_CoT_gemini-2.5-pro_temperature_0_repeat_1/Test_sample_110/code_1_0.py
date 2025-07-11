import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of a concrete encased composite column
    based on Eurocode 4 principles.
    """
    # 1. Define Material and Section Properties (in N, mm)
    print("--- 1. Material and Section Properties ---")
    # Storey and Column Geometry
    L = 4000.0  # Storey height in mm
    L_cr_factor = 1.0 # Effective length factor for pinned-pinned ends
    L_cr = L * L_cr_factor
    h_col = 400.0
    b_col = 400.0
    cover = 30.0

    # Structural Steel (S355, UC 254x254x132)
    f_y = 345.0  # Yield strength (N/mm^2)
    A_a = 168.0 * 100  # Area (mm^2)
    I_ay = 22500.0 * 1e4  # Second moment of area, y-axis (mm^4)
    I_az = 7530.0 * 1e4   # Second moment of area, z-axis (mm^4)
    E_a = 210000.0  # Modulus of Elasticity (N/mm^2)
    gamma_a = 1.0

    # Concrete (C40/50)
    f_ck = 40.0  # Characteristic cylinder strength (N/mm^2)
    gamma_c = 1.5
    alpha_cc = 0.85 # Factor for long term effects on compressive strength
    E_cm = 22 * ((f_ck + 8) / 10)**0.3 * 1000 # Secant modulus of elasticity (N/mm^2)

    # Reinforcing Steel (S500, 4xT16)
    f_sk = 500.0  # Characteristic yield strength (N/mm^2)
    rebar_dia = 16.0
    num_rebars = 4
    E_s = 200000.0  # Modulus of Elasticity (N/mm^2)
    gamma_s = 1.15
    
    print(f"Storey Height L = {L} mm")
    print(f"Steel (S355): f_y = {f_y} MPa, A_a = {A_a} mm^2, I_ay = {I_ay:.2e} mm^4, I_az = {I_az:.2e} mm^4")
    print(f"Concrete (C40/50): f_ck = {f_ck} MPa, E_cm = {E_cm:.0f} MPa")
    print(f"Rebar (S500): f_sk = {f_sk} MPa, Diameter = {rebar_dia} mm, Count = {num_rebars}")
    print("-" * 40 + "\n")

    # 2. Calculate Areas and Section Properties
    print("--- 2. Cross-Sectional Properties ---")
    A_c_gross = h_col * b_col
    A_s_one_bar = math.pi * (rebar_dia / 2)**2
    A_s = num_rebars * A_s_one_bar
    A_c = A_c_gross - A_a - A_s
    
    # Rebar position from centroid
    d_rebar = (h_col / 2) - cover - (rebar_dia / 2)
    # Second moment of area for rebars (using Parallel Axis Theorem, I_local is negligible)
    # All 4 bars contribute to both I_sy and I_sz as they are at the corners
    I_sy = num_rebars * A_s_one_bar * d_rebar**2
    I_sz = num_rebars * A_s_one_bar * d_rebar**2

    # Second moment of area for gross concrete section
    I_c_gross = (b_col * h_col**3) / 12
    
    print(f"Total Rebar Area A_s = {A_s:.2f} mm^2")
    print(f"Net Concrete Area A_c = {A_c:.2f} mm^2")
    print(f"Rebar Second Moment of Area I_sy = I_sz = {I_sy:.2e} mm^4")
    print(f"Gross Concrete Second Moment of Area I_c = {I_c_gross:.2e} mm^4")
    print("-" * 40 + "\n")

    # 3. Calculate Plastic Resistance (N_pl,Rd)
    print("--- 3. Plastic Resistance of Cross-Section (N_pl,Rd) ---")
    f_yd = f_y / gamma_a
    f_cd = alpha_cc * f_ck / gamma_c
    f_sd = f_sk / gamma_s
    alpha_c_pl = 1.0 # For encased sections as per EC4
    
    N_pl_Rd = (A_a * f_yd) + (alpha_c_pl * A_c * f_cd) + (A_s * f_sd)
    print(f"Design Strengths: f_yd={f_yd:.2f} MPa, f_cd={f_cd:.2f} MPa, f_sd={f_sd:.2f} MPa")
    print(f"N_pl,Rd = ({A_a:.0f}*{f_yd:.2f}) + ({alpha_c_pl:.1f}*{A_c:.0f}*{f_cd:.2f}) + ({A_s:.0f}*{f_sd:.2f})")
    print(f"N_pl,Rd = {N_pl_Rd/1000:.1f} kN")
    print("-" * 40 + "\n")

    # 4. Calculate Effective Flexural Stiffness (EI_eff)
    print("--- 4. Effective Flexural Stiffness (EI_eff) ---")
    K_e = 0.9 # Correction factor for steel contribution
    K_c = 0.6 # Factor for concrete contribution (short-term loading)

    EI_eff_y = K_e * E_a * I_ay + E_s * I_sy + K_c * E_cm * I_c_gross
    EI_eff_z = K_e * E_a * I_az + E_s * I_sz + K_c * E_cm * I_c_gross
    
    print(f"EI_eff,y = {EI_eff_y:.2e} Nmm^2")
    print(f"EI_eff,z = {EI_eff_z:.2e} Nmm^2")
    
    # 5. Determine Critical Buckling Axis
    if EI_eff_y < EI_eff_z:
        EI_eff = EI_eff_y
        axis = "y-y"
    else:
        EI_eff = EI_eff_z
        axis = "z-z"
    print(f"Critical buckling axis is {axis} with EI_eff = {EI_eff:.2e} Nmm^2")
    print("-" * 40 + "\n")

    # 6. Calculate Elastic Critical Buckling Load (N_cr)
    print("--- 5. Elastic Critical Buckling Load (N_cr) ---")
    N_cr = (math.pi**2 * EI_eff) / (L_cr**2)
    print(f"N_cr = (pi^2 * {EI_eff:.2e}) / {L_cr:.0f}^2")
    print(f"N_cr = {N_cr/1000:.1f} kN")
    print("-" * 40 + "\n")

    # 7. Calculate Relative Slenderness (lambda_bar)
    print("--- 6. Relative Slenderness (lambda_bar) ---")
    # Characteristic plastic resistance (gamma factors = 1.0, alpha_cc = 1.0)
    N_pl_Rk = (A_a * f_y) + (A_c * f_ck) + (A_s * f_sk)
    lambda_bar = math.sqrt(N_pl_Rk / N_cr)
    print(f"Characteristic Plastic Resistance N_pl,Rk = {N_pl_Rk/1000:.1f} kN")
    print(f"lambda_bar = sqrt({N_pl_Rk/1000:.1f} / {N_cr/1000:.1f}) = {lambda_bar:.3f}")
    print("-" * 40 + "\n")

    # 8. Determine Buckling Reduction Factor (chi)
    print("--- 7. Buckling Reduction Factor (chi) ---")
    # For encased I-section, weak axis (z-z) buckling, use buckling curve 'c'
    alpha = 0.49 # Imperfection factor for buckling curve 'c'
    phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    print(f"Using buckling curve 'c', imperfection factor alpha = {alpha}")
    print(f"phi = 0.5 * [1 + {alpha} * ({lambda_bar:.3f} - 0.2) + {lambda_bar:.3f}^2] = {phi:.3f}")
    print(f"chi = 1 / ({phi:.3f} + sqrt({phi:.3f}^2 - {lambda_bar:.3f}^2)) = {chi:.3f}")
    print("-" * 40 + "\n")

    # 9. Calculate Buckling Resistance (N_b,Rd)
    print("--- 8. Final Buckling Resistance (N_b,Rd) ---")
    N_b_Rd = chi * N_pl_Rd
    print(f"The buckling resistance is calculated as: N_b,Rd = chi * N_pl,Rd")
    print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd/1000:.1f} kN")
    print(f"Final Buckling Resistance N_b,Rd = {N_b_Rd/1000:.1f} kN")
    
    return N_b_Rd / 1000

# Execute the calculation and print the final answer
final_answer = calculate_buckling_resistance()
print(f"\n<<<The final buckling resistance is {final_answer:.1f} kN>>>")
<<<7630.0>>>