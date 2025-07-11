import math

def solve_buckling_resistance():
    """
    Calculates the buckling resistance of the concrete encased composite column.
    All calculations are in N and mm.
    """
    # Step 1: Define Material and Section Properties
    print("--- Step 1: Material and Section Properties ---")
    
    # Concrete C40/50
    f_ck = 40.0  # N/mm^2, characteristic cylinder strength
    gamma_c = 1.5
    f_cd = f_ck / gamma_c
    # Secant modulus of elasticity of concrete (EC2, Table 3.1)
    E_cm = 22 * ((f_ck + 8) / 10)**0.3 # GPa -> N/mm^2 is kN/mm^2 * 1000
    E_cm = E_cm * 1000 # N/mm^2
    print(f"Concrete C40/50: f_ck = {f_ck} MPa, f_cd = {f_cd:.2f} MPa, E_cm = {E_cm:.0f} MPa")

    # Structural Steel S355 (UC 254 x 254 x 132)
    f_y = 345.0  # N/mm^2
    gamma_a = 1.0 # Partial factor for structural steel
    f_yd = f_y / gamma_a
    A_a = 168.0 * 100  # cm^2 to mm^2
    I_ay = 22500.0 * 1e4  # cm^4 to mm^4
    I_az = 7530.0 * 1e4   # cm^4 to mm^4
    E_a = 210000.0 # N/mm^2
    print(f"Structural Steel S355: f_y = {f_y} MPa, A_a = {A_a} mm^2, I_ay = {I_ay:.1e} mm^4, I_az = {I_az:.1e} mm^4, E_a = {E_a} MPa")

    # Reinforcing Steel S500 (4xT16)
    f_sk = 500.0 # N/mm^2
    gamma_s = 1.15
    f_sd = f_sk / gamma_s
    d_rebar = 16.0 # mm
    n_rebars = 4
    A_s_one = math.pi * (d_rebar / 2)**2
    A_s = n_rebars * A_s_one
    E_s = 200000.0 # N/mm^2
    print(f"Reinforcing Steel S500: f_sk = {f_sk} MPa, A_s = {A_s:.2f} mm^2, E_s = {E_s} MPa")
    
    # Column Geometry
    h_c = 400.0 # mm
    b_c = 400.0 # mm
    cover = 30.0 # mm
    L = 4000.0 # mm, storey height
    # Buckling length factor k=1.0 (pinned-pinned) assumed as typical.
    k = 1.0
    L_cr = k * L
    print(f"Column Geometry: {h_c}x{b_c} mm, Length L_cr = {L_cr} mm")
    print("-" * 50)
    
    # Step 2: Calculate Areas and Inertias of Components
    print("--- Step 2: Component Areas and Inertias ---")
    A_c_gross = h_c * b_c
    A_c = A_c_gross - A_a - A_s # Net concrete area
    I_c_gross = (b_c * h_c**3) / 12
    
    # Rebar inertia
    rebar_dist = h_c / 2 - cover - d_rebar / 2
    I_sy = n_rebars * A_s_one * rebar_dist**2
    I_sz = n_rebars * A_s_one * rebar_dist**2
    print(f"Net Concrete Area A_c = {A_c:.2f} mm^2")
    print(f"Gross Concrete Inertia I_c = {I_c_gross:.2e} mm^4")
    print(f"Rebar Inertia I_sy = I_sz = {I_sy:.2e} mm^4")
    print("-" * 50)
    
    # Step 3: Calculate Plastic Resistance of Cross-Section (N_pl,Rd)
    print("--- Step 3: Plastic Resistance N_pl,Rd ---")
    # According to EC4 6.7.3.2(1), alpha_c = 1.0 for encased sections
    alpha_c = 1.0
    N_pl_Rd_a = A_a * f_yd
    N_pl_Rd_c = alpha_c * A_c * f_cd
    N_pl_Rd_s = A_s * f_sd
    N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s
    print(f"N_pl,Rd = A_a*f_yd + A_c*f_cd + A_s*f_sd")
    print(f"N_pl,Rd = {A_a:.0f}*{f_yd:.0f} + {A_c:.0f}*{f_cd:.2f} + {A_s:.0f}*{f_sd:.2f}")
    print(f"N_pl,Rd = {N_pl_Rd_a/1000:.1f} kN + {N_pl_Rd_c/1000:.1f} kN + {N_pl_Rd_s/1000:.1f} kN = {N_pl_Rd/1000:.1f} kN")
    print("-" * 50)

    # Step 4: Calculate Effective Flexural Stiffness ((EI)_eff)
    print("--- Step 4: Effective Flexural Stiffness (EI)_eff ---")
    K_e = 0.6 # Correction factor for concrete stiffness
    EI_eff_y = E_a * I_ay + E_s * I_sy + K_e * E_cm * I_c_gross
    EI_eff_z = E_a * I_az + E_s * I_sz + K_e * E_cm * I_c_gross
    print(f"(EI)_eff,y = E_a*I_ay + E_s*I_sy + K_e*E_cm*I_c = {EI_eff_y:.2e} N.mm^2")
    print(f"(EI)_eff,z = E_a*I_az + E_s*I_sz + K_e*E_cm*I_c = {EI_eff_z:.2e} N.mm^2")
    
    # Buckling occurs about the weaker axis (lower EI_eff)
    EI_eff = min(EI_eff_y, EI_eff_z)
    print(f"Buckling is critical about the z-z axis. (EI)_eff = {EI_eff:.2e} N.mm^2")
    print("-" * 50)
    
    # Step 5: Calculate Elastic Critical Buckling Load (N_cr)
    print("--- Step 5: Critical Buckling Load N_cr ---")
    N_cr = (math.pi**2 * EI_eff) / (L_cr**2)
    print(f"N_cr = pi^2 * (EI)_eff / L_cr^2")
    print(f"N_cr = pi^2 * {EI_eff:.2e} / {L_cr:.0f}^2 = {N_cr/1000:.1f} kN")
    print("-" * 50)

    # Step 6: Calculate Non-dimensional Slenderness (lambda_bar)
    print("--- Step 6: Non-dimensional Slenderness lambda_bar ---")
    N_pl_Rk = A_a * f_y + A_c * f_ck + A_s * f_sk
    print(f"N_pl,Rk = A_a*f_y + A_c*f_ck + A_s*f_sk = {N_pl_Rk/1000:.1f} kN")
    
    if N_cr > 0:
        lambda_bar = math.sqrt(N_pl_Rk / N_cr)
        print(f"lambda_bar = sqrt(N_pl,Rk / N_cr) = sqrt({N_pl_Rk/1000:.1f} / {N_cr/1000:.1f}) = {lambda_bar:.4f}")
    else:
        lambda_bar = 0
        print("N_cr is zero, cannot calculate slenderness.")
    print("-" * 50)
    
    # Step 7: Calculate Buckling Reduction Factor (chi)
    print("--- Step 7: Buckling Reduction Factor chi ---")
    # For encased I-section, buckling about z-z axis -> Buckling Curve 'c' (EC4 Table 6.5)
    # For curve 'c', imperfection factor alpha = 0.49 (EC3 Table 6.2)
    alpha = 0.49
    print(f"Using buckling curve 'c', imperfection factor alpha = {alpha}")
    
    phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
    print(f"phi = 0.5 * [1 + alpha*(lambda_bar-0.2) + lambda_bar^2]")
    print(f"phi = 0.5 * [1 + {alpha}*({lambda_bar:.4f}-0.2) + {lambda_bar:.4f}^2] = {phi:.4f}")

    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    if chi > 1.0:
      chi = 1.0
    print(f"chi = 1 / [phi + sqrt(phi^2 - lambda_bar^2)]")
    print(f"chi = 1 / [{phi:.4f} + sqrt({phi:.4f}^2 - {lambda_bar:.4f}^2)] = {chi:.4f}")
    print("-" * 50)
    
    # Step 8: Calculate Buckling Resistance (N_b,Rd)
    print("--- Step 8: Final Buckling Resistance N_b,Rd ---")
    N_b_Rd = chi * N_pl_Rd
    print(f"N_b,Rd = chi * N_pl,Rd = {chi:.4f} * {N_pl_Rd/1000:.1f} kN")
    print(f"The buckling resistance of the column is {N_b_Rd/1000:.1f} kN.")
    
    return N_b_Rd / 1000

if __name__ == "__main__":
    final_answer = solve_buckling_resistance()
    print(f"\n<<< {final_answer:.1f} >>>")
