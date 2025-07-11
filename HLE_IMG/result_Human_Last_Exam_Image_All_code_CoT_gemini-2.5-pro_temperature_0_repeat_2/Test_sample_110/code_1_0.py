import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of a concrete encased composite column
    based on Eurocode 4 (EN 1994-1-1).
    """
    # --- Step 1: Define Given Properties ---
    # Column properties
    L = 4000.0  # Storey height in mm
    L_cr_factor = 1.0 # Effective length factor for pinned-pinned ends
    L_cr = L * L_cr_factor # Effective buckling length in mm

    # Concrete properties (C40/50)
    f_ck = 40.0  # Characteristic cylinder strength in N/mm^2
    gamma_c = 1.5  # Partial safety factor for concrete
    alpha_cc = 0.85 # Factor for compressive strength
    K_c = 0.6 # Factor for effective stiffness of concrete (short-term loading)

    # Structural Steel properties (UC 254x254x132, S355)
    f_y = 345.0  # Yield strength in N/mm^2
    A_a_cm2 = 168.0  # Area in cm^2
    I_az_cm4 = 7530.0 # Second moment of area about z-z axis in cm^4
    E_a = 210000.0 # Modulus of Elasticity in N/mm^2
    gamma_a = 1.0 # Partial safety factor for steel
    K_e = 0.9 # Factor for effective stiffness of structural steel

    # Reinforcing Steel properties (4xT16, S500)
    num_rebars = 4
    d_rebar = 16.0 # Rebar diameter in mm
    f_sk = 500.0 # Characteristic yield strength in N/mm^2
    E_s = 200000.0 # Modulus of Elasticity in N/mm^2
    gamma_s = 1.15 # Partial safety factor for rebar

    # Section dimensions
    h_c = 400.0 # Concrete section height in mm
    b_c = 400.0 # Concrete section width in mm
    cover = 30.0 # Concrete cover in mm

    # --- Step 2: Convert Units and Calculate Derived Properties ---
    # Convert steel section properties to mm
    A_a = A_a_cm2 * 100  # mm^2
    I_az = I_az_cm4 * 10000 # mm^4

    # Design strengths
    f_yd = f_y / gamma_a
    f_cd = f_ck / gamma_c
    f_sd = f_sk / gamma_s

    # Concrete modulus of elasticity
    E_cm = 22 * ((f_ck + 8) / 10)**0.3

    # Rebar properties
    A_s1 = math.pi * (d_rebar / 2)**2 # Area of one rebar in mm^2
    A_s = num_rebars * A_s1 # Total area of rebars in mm^2

    # Concrete area
    A_c_gross = h_c * b_c
    A_c = A_c_gross - A_a - A_s

    # --- Step 3: Calculate Plastic Resistance to Compression (N_pl,Rd) ---
    N_pl_Rd_a = A_a * f_yd
    N_pl_Rd_c = alpha_cc * A_c * f_cd
    N_pl_Rd_s = A_s * f_sd
    N_pl_Rd = N_pl_Rd_a + N_pl_Rd_c + N_pl_Rd_s

    # --- Step 4: Calculate Effective Flexural Stiffness ((EI)_eff) for z-z axis ---
    d_rebar_centroid = h_c / 2 - cover - d_rebar / 2
    I_sz = num_rebars * A_s1 * d_rebar_centroid**2
    I_c_gross_z = (b_c * h_c**3) / 12
    I_cz = I_c_gross_z - I_az - I_sz
    EI_eff_z = K_e * E_a * I_az + E_s * I_sz + K_c * E_cm * I_cz

    # --- Step 5: Calculate Critical Buckling Load (N_cr) ---
    N_cr = (math.pi**2 * EI_eff_z) / (L_cr**2)

    # --- Step 6: Calculate Relative Slenderness (lambda_bar) ---
    N_pl_Rk = A_a * f_y + alpha_cc * A_c * f_ck + A_s * f_sk
    lambda_bar = math.sqrt(N_pl_Rk / N_cr)

    # --- Step 7: Determine Buckling Curve and Imperfection Factor (alpha) ---
    alpha = 0.49 # Buckling curve 'c' for encased I-section about minor axis

    # --- Step 8: Calculate Reduction Factor (chi) ---
    phi = 0.5 * (1 + alpha * (lambda_bar - 0.2) + lambda_bar**2)
    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    chi = min(chi, 1.0)

    # --- Step 9: Calculate Buckling Resistance (N_b,Rd) ---
    N_b_Rd = chi * N_pl_Rd

    # --- Step 10: Print the results ---
    print("--- Calculation of Buckling Resistance of Composite Column ---")
    print("\nStep 1: Plastic Resistance of the Cross-Section (N_pl,Rd)")
    print(f"N_pl,Rd = (A_a*f_yd) + (0.85*A_c*f_cd) + (A_s*f_sd)")
    print(f"N_pl,Rd = ({A_a:.1f} * {f_yd:.1f}) + (0.85 * {A_c:.1f} * {f_cd:.2f}) + ({A_s:.1f} * {f_sd:.2f})")
    print(f"N_pl,Rd = {N_pl_Rd_a/1000:.1f} kN + {N_pl_Rd_c/1000:.1f} kN + {N_pl_Rd_s/1000:.1f} kN = {N_pl_Rd/1000:.1f} kN")

    print("\nStep 2: Critical Buckling Load (N_cr)")
    print(f"Effective Flexural Stiffness (EI)_eff = {EI_eff_z:.3e} Nmm^2")
    print(f"N_cr = pi^2 * (EI)_eff / L_cr^2 = {N_cr/1000:.1f} kN")

    print("\nStep 3: Relative Slenderness (lambda_bar)")
    print(f"Characteristic Plastic Resistance N_pl,Rk = {N_pl_Rk/1000:.1f} kN")
    print(f"lambda_bar = sqrt(N_pl,Rk / N_cr) = sqrt({N_pl_Rk/1000:.1f} / {N_cr/1000:.1f}) = {lambda_bar:.4f}")

    print("\nStep 4: Buckling Reduction Factor (chi)")
    print(f"For buckling curve 'c', imperfection factor alpha = {alpha}")
    print(f"chi = {chi:.4f}")

    print("\nStep 5: Final Buckling Resistance (N_b,Rd)")
    print(f"The buckling resistance is calculated as: N_b,Rd = chi * N_pl,Rd")
    print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd/1000:.1f}")
    print(f"N_b,Rd = {N_b_Rd/1000:.1f} kN")

    # Return the final answer for the platform
    return N_b_Rd/1000

if __name__ == '__main__':
    final_answer = calculate_buckling_resistance()
    # print(f"\n<<< {final_answer:.1f} >>>") # For internal check
    
calculate_buckling_resistance()