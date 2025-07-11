import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of the concrete encased composite column.
    All calculations are in N and mm.
    """
    # 1. Define Material and Geometric Properties
    # Storey height
    L = 4000.0  # mm
    
    # Concrete C40/50
    f_ck = 40.0  # N/mm^2
    gamma_c = 1.5
    
    # Structural Steel UC 254x254x132 S355
    f_y = 345.0  # N/mm^2
    gamma_a = 1.0
    A_a = 168.0 * 100  # mm^2 (168 cm^2)
    I_ay = 22500.0 * 1e4 # mm^4 (strong axis)
    I_az = 7530.0 * 1e4  # mm^4 (weak axis)
    
    # Reinforcement S500
    f_sk = 500.0  # N/mm^2
    gamma_s = 1.15
    num_rebars = 4
    rebar_dia = 16.0  # mm
    
    # Section geometry
    h_col = 400.0 # mm
    b_col = 400.0 # mm
    cover = 30.0 # mm
    
    # --- Step 2: Calculate Material Design Strengths & Moduli ---
    f_yd = f_y / gamma_a
    # According to EN 1994-1-1 6.7.3.2(3), Npl,Rd = Aa*fyd + 0.85*Ac*fcd + As*fsd
    # where fcd = fck / gamma_c. We use alpha_cc=1.0 as long term effects are ignored.
    f_cd = f_ck / gamma_c
    f_sd = f_sk / gamma_s
    
    # Elastic Moduli
    E_a = 210000.0  # N/mm^2
    E_s = 210000.0  # N/mm^2
    E_cm = 22 * ((f_ck + 8) / 10)**0.3 * 1000 # N/mm^2 (Secant modulus from EC2)

    # --- Step 3: Calculate Cross-Sectional Areas and Plastic Resistance N_pl,Rd ---
    A_rebar = math.pi * (rebar_dia / 2)**2
    A_s = num_rebars * A_rebar
    A_c_gross = h_col * b_col
    A_c = A_c_gross - A_a - A_s
    
    # Plastic compression resistance N_pl,Rd
    N_pl_Rd = (A_a * f_yd) + (0.85 * A_c * f_cd) + (A_s * f_sd)

    # --- Step 4: Calculate Effective Flexural Stiffness (EI)_eff ---
    # Second moment of area of concrete (uncracked gross section)
    I_c = b_col * h_col**3 / 12
    
    # Second moment of area of reinforcement
    rebar_dist_from_center = (h_col / 2) - cover - (rebar_dia / 2)
    # I_s = sum(A_si * d_i^2), where d is distance from centroidal axis.
    # For both y and z axes, all 4 bars are at distance `rebar_dist_from_center`.
    I_s = num_rebars * A_rebar * rebar_dist_from_center**2

    # Stiffness Correction Factor (creep ignored)
    K_e = 0.6
    
    # Effective flexural stiffness (EI)_eff for weak (z-z) axis, which governs
    EI_eff_z = E_a * I_az + E_s * I_s + K_e * E_cm * I_c

    # --- Step 5: Calculate Critical Buckling Load (N_cr) ---
    # Assume pinned-pinned column, Buckling length L_cr = L
    L_cr = L
    N_cr = (math.pi**2 * EI_eff_z) / L_cr**2
    
    # --- Step 6: Calculate Relative Slenderness (lambda_bar) ---
    # Characteristic plastic resistance N_pl,Rk (using characteristic strengths)
    N_pl_Rk = (A_a * f_y) + (0.85 * A_c * f_ck) + (A_s * f_sk)
    
    # Relative slenderness
    if N_cr > 0:
      lambda_bar = math.sqrt(N_pl_Rk / N_cr)
    else:
      lambda_bar = 0

    # --- Step 7: Determine Buckling Reduction Factor (chi) ---
    # For weak axis (z-z) buckling of encased I-section, use buckling curve 'b'.
    alpha = 0.34 # Imperfection factor for buckling curve 'b'
    
    phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    
    # --- Step 8: Calculate Buckling Resistance (N_b,Rd) ---
    N_b_Rd = chi * N_pl_Rd
    
    # --- Step 9: Final Output ---
    N_pl_Rd_kN = N_pl_Rd / 1000
    N_b_Rd_kN = N_b_Rd / 1000
    
    print("--- Calculation Summary ---")
    print(f"Plastic Compression Resistance (N_pl,Rd): {N_pl_Rd_kN:.2f} kN")
    print(f"Effective Flexural Stiffness (EI_eff,z): {EI_eff_z / 1e12:.2f} x 10^12 N.mm^2")
    print(f"Critical Buckling Load (N_cr): {N_cr / 1000:.2f} kN")
    print(f"Relative Slenderness (λ_bar): {lambda_bar:.4f}")
    print(f"Buckling Reduction Factor (χ): {chi:.4f}")
    print("\n--- Final Calculation ---")
    print(f"Buckling Resistance N_b,Rd = χ * N_pl,Rd")
    print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd_kN:.2f} kN")
    print(f"N_b,Rd = {N_b_Rd_kN:.2f} kN")

    return N_b_Rd_kN

if __name__ == '__main__':
    buckling_resistance = calculate_buckling_resistance()
    # The final answer in the required format
    # The calculated value is 8189.60, which rounds to 8190
    print(f"<<<{buckling_resistance:.1f}>>>")
