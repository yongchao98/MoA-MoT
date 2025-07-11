import math

def solve_buckling_resistance():
    """
    Calculates the buckling resistance of a concrete encased composite column
    according to Eurocode 4 (EN 1994-1-1).
    """
    # --- 1. Given Data & Material Properties ---
    # Geometry
    storey_height = 4000.0  # mm
    h_total = 400.0  # mm
    b_total = 400.0  # mm
    cover = 30.0  # mm
    rebar_dia = 16.0  # mm
    num_rebars = 4

    # Concrete C40/50
    f_ck = 40.0  # MPa (N/mm^2)
    gamma_c = 1.5
    E_cm = 35000.0  # MPa

    # Structural Steel UC 254x254x132 S355
    f_y = 345.0  # MPa
    A_a = 168.0 * 100  # mm^2 from cm^2
    I_ay = 22500.0 * 1e4  # mm^4 from cm^4 (Strong axis)
    I_az = 7530.0 * 1e4   # mm^4 from cm^4 (Weak axis)
    gamma_a = 1.0
    E_a = 210000.0  # MPa

    # Reinforcing Steel S500
    f_sk = 500.0  # MPa
    gamma_s = 1.15
    E_s = 200000.0  # MPa

    print("--- 1. Initial Properties ---")
    print(f"Storey height, L = {storey_height} mm")
    print(f"Concrete: C40/50, f_ck = {f_ck} MPa")
    print(f"Structural Steel: S355, f_y = {f_y} MPa")
    print(f"Reinforcing Steel: S500, f_sk = {f_sk} MPa\n")

    # --- 2. Sectional Properties Calculations ---
    print("--- 2. Sectional Properties ---")
    # Gross Concrete Section
    A_gross = h_total * b_total
    I_c_gross = b_total * h_total**3 / 12  # Same for both axes
    
    # Reinforcement
    A_s1 = math.pi * (rebar_dia**2) / 4
    A_s = num_rebars * A_s1
    dist_rebar_center = h_total/2 - cover
    I_s = num_rebars * A_s1 * dist_rebar_center**2

    # Net Concrete Area
    A_c_net = A_gross - A_a - A_s

    print(f"Area of Steel Section, A_a = {A_a:.2f} mm^2")
    print(f"Area of Reinforcement, A_s = {A_s:.2f} mm^2")
    print(f"Net Area of Concrete, A_c,net = {A_c_net:.2f} mm^2")
    print(f"Inertia of Gross Concrete Section, I_c = {I_c_gross:.2e} mm^4")
    print(f"Inertia of Reinforcement (about centroid), I_s = {I_s:.2e} mm^4\n")

    # --- 3. Plastic Resistance (N_pl,Rd) ---
    print("--- 3. Plastic Resistance Calculation (N_pl,Rd) ---")
    f_yd = f_y / gamma_a
    f_sd = f_sk / gamma_s
    f_cd = f_ck / gamma_c

    N_a_pl_Rd = A_a * f_yd
    N_s_pl_Rd = A_s * f_sd
    N_c_pl_Rd = 0.85 * A_c_net * f_cd

    N_pl_Rd = N_a_pl_Rd + N_s_pl_Rd + N_c_pl_Rd
    
    print(f"Design resistance from steel section = {N_a_pl_Rd/1000:.2f} kN")
    print(f"Design resistance from rebar = {N_s_pl_Rd/1000:.2f} kN")
    print(f"Design resistance from concrete = {N_c_pl_Rd/1000:.2f} kN")
    print(f"Total Plastic Resistance, N_pl,Rd = {N_pl_Rd/1000:.2f} kN\n")

    # --- 4. Effective Flexural Stiffness (EI_eff) ---
    print("--- 4. Effective Flexural Stiffness Calculation (EI_eff) ---")
    # Per EN1994-1-1 Eq 6.40, K_e = 0.6 is a calibration factor. This is used
    # even when ignoring long-term effects (creep), as it also accounts
    # for model uncertainties. Ignoring long-term effects means we use E_cm.
    K_e = 0.6
    
    # Stiffness about y-y axis (strong axis of steel section)
    EI_eff_y = E_a * I_ay + K_e * E_cm * I_c_gross + E_s * I_s
    
    # Stiffness about z-z axis (weak axis of steel section)
    EI_eff_z = E_a * I_az + K_e * E_cm * I_c_gross + E_s * I_s

    EI_eff = min(EI_eff_y, EI_eff_z)
    
    print(f"EI_eff,y = {EI_eff_y:.2e} Nmm^2")
    print(f"EI_eff,z = {EI_eff_z:.2e} Nmm^2")
    print(f"Governing Effective Stiffness, EI_eff = {EI_eff:.2e} Nmm^2 (buckling about weak axis)\n")

    # --- 5. Critical Buckling Load (N_cr) ---
    print("--- 5. Critical Buckling Load Calculation (N_cr) ---")
    # Assume pinned-pinned column, buckling length L_cr = 1.0 * L
    L_cr = 1.0 * storey_height
    N_cr = (math.pi**2 * EI_eff) / (L_cr**2)
    
    print(f"Buckling Length, L_cr = {L_cr:.0f} mm")
    print(f"Elastic Critical Buckling Load, N_cr = {N_cr/1000:.2f} kN\n")
    
    # --- 6. Relative Slenderness (lambda_bar) ---
    print("--- 6. Relative Slenderness Calculation (lambda_bar) ---")
    # Need characteristic plastic resistance N_pl,Rk
    N_a_pl_Rk = A_a * f_y
    N_s_pl_Rk = A_s * f_sk
    N_c_pl_Rk = 0.85 * A_c_net * f_ck
    N_pl_Rk = N_a_pl_Rk + N_s_pl_Rk + N_c_pl_Rk

    lambda_bar = math.sqrt(N_pl_Rk / N_cr)
    
    print(f"Characteristic Plastic Resistance, N_pl,Rk = {N_pl_Rk/1000:.2f} kN")
    print(f"Relative Slenderness, lambda_bar = {lambda_bar:.4f}\n")

    # --- 7. Buckling Reduction Factor (chi) ---
    print("--- 7. Buckling Reduction Factor Calculation (chi) ---")
    # For concrete encased I-section, weak axis buckling -> buckling curve 'b'
    # From EN1993-1-1 Table 6.1, for curve 'b', alpha = 0.34
    alpha = 0.34
    
    phi_term = alpha * (lambda_bar - 0.2) if lambda_bar > 0.2 else 0
    phi = 0.5 * (1 + phi_term + lambda_bar**2)
    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    chi = min(chi, 1.0)
    
    print(f"Buckling curve 'b' used, imperfection factor, alpha = {alpha}")
    print(f"Buckling Reduction Factor, chi = {chi:.4f}\n")

    # --- 8. Final Buckling Resistance (N_b,Rd) ---
    print("--- 8. Final Buckling Resistance Calculation (N_b,Rd) ---")
    N_b_Rd = chi * N_pl_Rd
    
    print(f"The design buckling resistance is calculated as:")
    print(f"N_b,Rd = chi * N_pl,Rd")
    print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd/1000:.2f} kN")
    
    final_answer = N_b_Rd / 1000
    print(f"Final Design Buckling Resistance, N_b,Rd = {final_answer:.0f} kN")
    
if __name__ == '__main__':
    solve_buckling_resistance()
    print("\n<<<8190>>>")
