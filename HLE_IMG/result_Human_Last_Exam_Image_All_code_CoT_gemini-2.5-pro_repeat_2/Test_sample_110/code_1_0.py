import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of the concrete encased composite column.
    The method follows the principles of EN 1994-1-1.
    """
    # --- 1. Define constants and input parameters ---
    print("--- Input Parameters ---")
    # Steel Section (S355)
    f_y = 345.0  # N/mm^2, Yield strength
    E_a = 210000.0  # N/mm^2, Modulus of Elasticity
    gamma_a = 1.0  # Partial safety factor
    A_a = 168.0 * 100  # cm^2 to mm^2, Area
    I_az = 7530.0 * 1e4 # cm^4 to mm^4, Second moment of area (weak axis)
    print(f"Steel (S355): f_y={f_y} N/mm^2, A_a={A_a} mm^2, I_az={I_az:.2e} mm^4")

    # Concrete (C40/50)
    f_ck = 40.0  # N/mm^2, Characteristic compressive strength
    gamma_c = 1.5  # Partial safety factor
    print(f"Concrete (C40/50): f_ck={f_ck} N/mm^2")

    # Reinforcement (S500)
    f_sk = 500.0  # N/mm^2, Characteristic yield strength
    E_s = 200000.0  # N/mm^2, Modulus of Elasticity
    gamma_s = 1.15 # Partial safety factor
    d_rebar = 16.0 # mm, Diameter
    n_rebar = 4    # Number of rebars
    print(f"Rebars (S500): 4xT{int(d_rebar)}, f_sk={f_sk} N/mm^2")

    # Section & Column Dimensions
    h_c = 400.0 # mm, Overall height
    b_c = 400.0 # mm, Overall width
    cover = 30.0 # mm, Concrete cover
    L = 4000.0 # mm, Storey height
    k_L = 1.0  # Effective length factor (assuming pinned ends)
    L_cr = k_L * L # mm, Buckling length
    print(f"Column: {h_c}x{b_c} mm, Length L_cr={L_cr} mm")
    print("-" * 25 + "\n")

    # --- 2. Calculate Derived Section & Material Properties ---
    print("--- Step 1: Calculate Section & Material Properties ---")
    # Concrete modulus of elasticity
    E_cm = 22000 * ((f_ck + 8) / 10)**0.3
    # Reinforcement area
    A_s1 = math.pi * (d_rebar / 2)**2
    A_s = n_rebar * A_s1
    # Net concrete area
    A_c = h_c * b_c - A_a - A_s
    # Second moment of area of reinforcement (about weak z-z axis)
    dist_rebar_z = b_c / 2 - cover - d_rebar / 2
    I_sz = n_rebar * A_s1 * dist_rebar_z**2
    # Second moment of area of gross concrete (about weak z-z axis)
    I_cz_gross = b_c * h_c**3 / 12
    print(f"Net Concrete Area A_c = {A_c:.2f} mm^2")
    print(f"Total Rebar Area A_s = {A_s:.2f} mm^2")
    print(f"Rebar Second Moment of Area I_sz = {I_sz:.2e} mm^4")
    print("-" * 25 + "\n")

    # --- 3. Calculate Plastic Resistance (N_pl,Rd) ---
    print("--- Step 2: Calculate Design Plastic Resistance (N_pl,Rd) ---")
    f_yd = f_y / gamma_a
    f_sd = f_sk / gamma_s
    f_cd = f_ck / gamma_c
    alpha_c = 1.0 # For encased sections

    N_pl_Rd_a = A_a * f_yd
    N_pl_Rd_s = A_s * f_sd
    N_pl_Rd_c = alpha_c * A_c * f_cd
    N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s
    print(f"N_pl,Rd = (A_a*f_y/γ_a) + (A_c*f_ck/γ_c) + (A_s*f_sk/γ_s)")
    print(f"N_pl,Rd = ({N_pl_Rd_a/1000:.1f}) + ({N_pl_Rd_c/1000:.1f}) + ({N_pl_Rd_s/1000:.1f})")
    print(f"Total N_pl,Rd = {N_pl_Rd / 1000:.2f} kN")
    print("-" * 25 + "\n")

    # --- 4. Calculate Effective Flexural Stiffness ((EI)_eff) ---
    print("--- Step 3: Calculate Effective Flexural Stiffness ((EI)_eff,z) ---")
    K_e = 0.9  # Correction factor
    K_c = 0.6  # Factor for concrete contribution
    
    EI_a = E_a * I_az
    EI_s = E_s * I_sz
    EI_c = K_c * E_cm * I_cz_gross
    EI_eff_z = K_e * (EI_a + EI_s + EI_c)
    print(f"(EI)_eff,z = K_e * (E_a*I_az + E_s*I_sz + K_c*E_cm*I_cz)")
    print(f"(EI)_eff,z = {K_e} * ({EI_a:.2e} + {EI_s:.2e} + {EI_c:.2e})")
    print(f"(EI)_eff,z = {EI_eff_z:.2e} Nmm^2")
    print("-" * 25 + "\n")

    # --- 5. Calculate Critical Buckling Load (N_cr) ---
    print("--- Step 4: Calculate Critical Buckling Load (N_cr) ---")
    N_cr = (math.pi**2 * EI_eff_z) / (L_cr**2)
    print(f"N_cr = (pi^2 * (EI)_eff,z) / L_cr^2")
    print(f"N_cr = (pi^2 * {EI_eff_z:.2e}) / {L_cr}^2 = {N_cr / 1000:.2f} kN")
    print("-" * 25 + "\n")

    # --- 6. Calculate Slenderness and Reduction Factor (chi) ---
    print("--- Step 5: Calculate Reduction Factor (chi) ---")
    # Characteristic plastic resistance (for slenderness calculation)
    N_pl_Rk = A_a * f_y + A_s * f_sk + A_c * f_ck
    # Relative slenderness
    lambda_bar = math.sqrt(N_pl_Rk / N_cr)
    print(f"Relative Slenderness λ_bar = sqrt(N_pl,Rk / N_cr) = {lambda_bar:.4f}")

    # Buckling curve 'c' -> imperfection factor alpha
    alpha_imp = 0.49
    # Intermediate factor phi
    phi = 0.5 * (1 + alpha_imp * (lambda_bar - 0.2) + lambda_bar**2)
    print(f"Factor phi (for buckling curve 'c', α=0.49) = {phi:.4f}")
    
    # Reduction factor chi
    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    chi = min(chi, 1.0) # Cannot be greater than 1.0
    print(f"Reduction Factor χ = 1 / (phi + sqrt(phi^2 - λ_bar^2)) = {chi:.4f}")
    print("-" * 25 + "\n")

    # --- 7. Calculate Buckling Resistance (N_b,Rd) ---
    print("--- Step 6: Calculate Final Buckling Resistance (N_b,Rd) ---")
    N_b_Rd = chi * N_pl_Rd
    print(f"N_b,Rd = χ * N_pl,Rd")
    print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd / 1000:.2f}")
    print(f"\nThe buckling resistance of the column is {N_b_Rd / 1000:.1f} kN.")
    
    return N_b_Rd / 1000

if __name__ == '__main__':
    final_answer = calculate_buckling_resistance()
