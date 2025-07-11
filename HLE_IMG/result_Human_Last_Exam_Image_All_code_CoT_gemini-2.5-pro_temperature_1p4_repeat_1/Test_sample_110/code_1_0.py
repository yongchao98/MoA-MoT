import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of a concrete encased composite column
    according to EN 1994-1-1.
    """
    # --- 1. Inputs and Material Properties ---
    # Partial safety factors
    gamma_a = 1.0  # For structural steel
    gamma_s = 1.15 # For reinforcing steel
    gamma_c = 1.5  # For concrete

    # Material characteristic strengths
    f_y = 345  # N/mm^2 (Yield strength of S355 steel)
    f_sk = 500 # N/mm^2 (Yield strength of S500 rebar)
    f_ck = 40  # N/mm^2 (Characteristic cylinder strength of C40/50 concrete)

    # Material design strengths
    f_yd = f_y / gamma_a
    f_sd = f_sk / gamma_s
    alpha_cc = 0.85 # Factor for long term effects on concrete strength
    f_cd = alpha_cc * f_ck / gamma_c

    # Moduli of Elasticity
    E_a = 210000  # N/mm^2 (Steel section)
    E_s = 200000  # N/mm^2 (Reinforcement)
    E_cm = 22000 * ((f_ck + 8) / 10)**0.3 # N/mm^2 (Secant modulus of concrete)

    # --- 2. Section Properties ---
    # Column dimensions
    h_col = 400  # mm
    b_col = 400  # mm
    cover = 30 # mm

    # Steel Section (UC 254 x 254 x 132)
    A_a = 16800  # mm^2 (Area)
    I_az = 7530 * 1e4 # mm^4 (Second moment of area about weaker z-z axis)

    # Reinforcement (4 x T16)
    d_bar = 16  # mm (diameter)
    A_s1 = math.pi * (d_bar / 2)**2 # Area of one bar
    A_s = 4 * A_s1
    # Distance from section centerline to rebar centerline
    d_rebar_z = h_col / 2 - cover - d_bar / 2
    I_sz = A_s * d_rebar_z**2 # Second moment of area for rebars about z-z axis

    # Concrete
    A_c_gross = h_col * b_col
    A_c = A_c_gross - A_a - A_s
    I_c = b_col * h_col**3 / 12

    print("--- Intermediate Calculation Values ---")
    print(f"Design Strength of Structural Steel (fyd): {f_yd:.2f} N/mm^2")
    print(f"Design Strength of Rebar (fsd): {f_sd:.2f} N/mm^2")
    print(f"Design Strength of Concrete (fcd): {f_cd:.2f} N/mm^2")
    print(f"Total Area of Reinforcement (As): {A_s:.2f} mm^2")
    print(f"Net Area of Concrete (Ac): {A_c:.2f} mm^2\n")

    # --- 3. Plastic Resistance of Cross-Section (N_pl,Rd) ---
    N_pl_Rd = A_a * f_yd + A_c * f_cd + A_s * f_sd
    print(f"Plastic Resistance (N_pl,Rd) = {A_a:.0f}*{f_yd:.2f} + {A_c:.2f}*{f_cd:.2f} + {A_s:.2f}*{f_sd:.2f}")
    print(f"N_pl,Rd = {N_pl_Rd/1000:.2f} kN\n")

    # --- 4. Effective Flexural Stiffness ((EI)_eff) ---
    K_0 = 0.9  # Calibration factor
    K_e = 0.6  # Factor for concrete stiffness (short-term loading)
    EI_eff = K_0 * (E_a * I_az + E_s * I_sz + K_e * E_cm * I_c)
    print(f"Effective Flexural Stiffness (EI_eff) = {K_0} * ({E_a:.0f}*{I_az:.2e} + {E_s:.0f}*{I_sz:.2e} + {K_e}*{E_cm:.2f}*{I_c:.2e})")
    print(f"(EI)_eff = {EI_eff:.2e} Nmm^2\n")

    # --- 5. Critical Buckling Load (N_cr) ---
    L = 4000  # mm (Storey height)
    L_cr = 1.0 * L # Buckling length (assuming pinned-pinned ends)
    N_cr = (math.pi**2 * EI_eff) / (L_cr**2)
    print(f"Critical Buckling Load (N_cr) = (pi^2 * {EI_eff:.2e}) / {L_cr:.0f}^2")
    print(f"N_cr = {N_cr/1000:.2f} kN\n")

    # --- 6. Relative Slenderness (lambda_bar) ---
    # Characteristic plastic resistance N_pl,Rk (using characteristic strengths)
    N_pl_Rk = A_a * f_y + 0.85 * A_c * f_ck + A_s * f_sk
    lambda_bar = math.sqrt(N_pl_Rk / N_cr)
    print(f"Characteristic Plastic Resistance (N_pl,Rk): {N_pl_Rk/1000:.2f} kN")
    print(f"Relative Slenderness (lambda_bar) = sqrt({N_pl_Rk:.2f} / {N_cr:.2f})")
    print(f"lambda_bar = {lambda_bar:.4f}\n")

    # --- 7. Buckling Reduction Factor (chi) ---
    # For UC section buckling about z-z axis, use buckling curve 'c'
    alpha = 0.49 # Imperfection factor for buckling curve 'c'
    Phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
    chi = 1 / (Phi + math.sqrt(Phi**2 - lambda_bar**2))
    
    # Ensure chi is not greater than 1
    chi = min(chi, 1.0)
    print(f"Imperfection Factor (alpha): {alpha}")
    print(f"Phi = 0.5 * (1 + {alpha} * ({lambda_bar:.4f} - 0.2) + {lambda_bar:.4f}^2)")
    print(f"Phi = {Phi:.4f}")
    print(f"Buckling Reduction Factor (chi) = 1 / ({Phi:.4f} + sqrt({Phi:.4f}^2 - {lambda_bar:.4f}^2))")
    print(f"chi = {chi:.4f}\n")
    
    # --- 8. Buckling Resistance (N_b,Rd) ---
    N_b_Rd = chi * N_pl_Rd
    
    print("--- Final Calculation ---")
    print(f"Buckling Resistance (N_b,Rd) = chi * N_pl,Rd")
    print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd/1000:.2f} kN")
    print(f"N_b,Rd = {N_b_Rd/1000:.2f} kN")
    
    return N_b_Rd / 1000

# Execute the calculation and print the final answer
buckling_resistance_kN = calculate_buckling_resistance()
print("\nFinal Answer:")
print(f"The buckling resistance of the column is {buckling_resistance_kN:.1f} kN.")
print(f"<<<{buckling_resistance_kN:.1f}>>>")
