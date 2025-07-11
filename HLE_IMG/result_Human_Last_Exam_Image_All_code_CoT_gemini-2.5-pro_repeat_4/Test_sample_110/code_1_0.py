import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of a concrete encased composite column
    based on Eurocode 4.
    """
    print("Step 1: Define Material and Sectional Properties\n")

    # --- Material Properties (with partial safety factors) ---
    # Structural Steel (UC 254x254x132 S355)
    f_y = 345.0      # N/mm^2 (Yield strength)
    E_a = 210000.0   # N/mm^2 (Modulus of Elasticity)
    gamma_a = 1.0    # Partial safety factor for steel

    # Concrete (C40/50)
    f_ck = 40.0      # N/mm^2 (Characteristic cylinder strength)
    gamma_c = 1.5    # Partial safety factor for concrete
    # Mean modulus of elasticity for concrete (EC2, Table 3.1)
    E_cm = 22 * ((f_ck + 8) / 10)**0.3 * 1000 # N/mm^2

    # Reinforcing Steel (T16 S500)
    d_rebar = 16.0   # mm (Diameter)
    f_sk = 500.0     # N/mm^2 (Characteristic yield strength)
    E_s = 200000.0   # N/mm^2 (Modulus of Elasticity)
    gamma_s = 1.15   # Partial safety factor for reinforcement

    # --- Sectional Properties ---
    # Steel Section (UC 254 x 254 x 132)
    A_a = 168.0 * 100  # mm^2 (from 168 cm^2)
    I_az = 7530.0 * 1e4 # mm^4 (from 7530 cm^4, weak axis)

    # Concrete Section
    h_c = 400.0      # mm
    b_c = 400.0      # mm
    A_c_gross = h_c * b_c
    I_cz = (h_c * b_c**3) / 12  # mm^4

    # Reinforcement Section
    n_rebars = 4
    A_s1 = math.pi * (d_rebar / 2)**2
    A_s = n_rebars * A_s1
    
    # Rebar position
    cover = 30.0 # mm
    d_from_center = h_c / 2 - cover - d_rebar / 2
    # Second moment of area of rebars about section's z-axis
    I_sz = A_s * d_from_center**2

    # Net Concrete Area
    A_c = A_c_gross - A_a - A_s
    
    print(f"Steel section area (A_a): {A_a:.2f} mm^2")
    print(f"Concrete net area (A_c): {A_c:.2f} mm^2")
    print(f"Rebar area (A_s): {A_s:.2f} mm^2\n")

    print("Step 2: Calculate Plastic Resistance (N_pl,Rd)\n")
    f_yd = f_y / gamma_a
    f_cd = f_ck / gamma_c
    f_sd = f_sk / gamma_s
    alpha_c = 0.85 # Factor for encased sections

    Npl_Rd_a = A_a * f_yd
    Npl_Rd_c = alpha_c * A_c * f_cd
    Npl_Rd_s = A_s * f_sd
    Npl_Rd = Npl_Rd_a + Npl_Rd_c + Npl_Rd_s

    print(f"Design Plastic Resistance N_pl,Rd = {Npl_Rd/1000:.2f} kN\n")

    print("Step 3: Calculate Effective Flexural Stiffness (EI_eff) about the weaker z-axis\n")
    K_e = 0.6  # Factor for short-term loading (creep ignored)
    EI_eff_z = E_a * I_az + K_e * E_cm * I_cz + E_s * I_sz
    
    print(f"Effective Flexural Stiffness (EI_eff,z): {EI_eff_z:.2e} Nmm^2\n")

    print("Step 4: Calculate Critical Buckling Load (N_cr)\n")
    storey_height = 4000.0  # mm
    L_cr = 1.0 * storey_height  # Effective length for pinned ends
    N_cr_z = (math.pi**2 * EI_eff_z) / L_cr**2

    print(f"Effective Length (L_cr): {L_cr:.0f} mm")
    print(f"Critical Buckling Load (N_cr,z): {N_cr_z/1000:.2f} kN\n")

    print("Step 5: Calculate Non-dimensional Slenderness (λ_bar)\n")
    # Characteristic plastic resistance (using characteristic strengths)
    Npl_Rk = A_a * f_y + alpha_c * A_c * f_ck + A_s * f_sk
    lambda_bar_z = math.sqrt(Npl_Rk / N_cr_z)
    
    print(f"Characteristic Plastic Resistance (N_pl,Rk): {Npl_Rk/1000:.2f} kN")
    print(f"Non-dimensional Slenderness (λ_bar_z): {lambda_bar_z:.4f}\n")

    print("Step 6: Determine Buckling Reduction Factor (χ)\n")
    # For encased I-section buckling about z-z axis, use buckling curve 'b'
    alpha = 0.34  # Imperfection factor for curve 'b'
    Phi_z = 0.5 * (1 + alpha * (lambda_bar_z - 0.2) + lambda_bar_z**2)
    chi_z = 1 / (Phi_z + math.sqrt(Phi_z**2 - lambda_bar_z**2))
    chi_z = min(chi_z, 1.0) # chi cannot be greater than 1.0

    print(f"Using Buckling Curve 'b' (alpha = {alpha})")
    print(f"Intermediate value Phi: {Phi_z:.4f}")
    print(f"Buckling Reduction Factor (χ_z): {chi_z:.4f}\n")

    print("Step 7: Calculate Buckling Resistance (N_b,Rd)\n")
    N_b_Rd = chi_z * Npl_Rd

    print("The final equation for Buckling Resistance is:")
    print("N_b,Rd = χ * N_pl,Rd")
    print("N_b,Rd = χ * (A_a * f_yd + α_c * A_c * f_cd + A_s * f_sd)")
    print(f"N_b,Rd = {chi_z:.4f} * ({A_a:.0f} * {f_yd:.2f} + {alpha_c} * {A_c:.2f} * {f_cd:.2f} + {A_s:.2f} * {f_sd:.2f})")
    print(f"N_b,Rd = {chi_z:.4f} * ({Npl_Rd_a:.0f} + {Npl_Rd_c:.0f} + {Npl_Rd_s:.0f})")
    print(f"N_b,Rd = {chi_z:.4f} * {Npl_Rd:.0f} N")
    print(f"N_b,Rd = {N_b_Rd/1000:.1f} kN\n")
    
    print("--- FINAL ANSWER ---")
    print(f"The buckling resistance of the column is {N_b_Rd/1000:.1f} kN.")
    
    return N_b_Rd/1000

# Execute the calculation
final_answer = calculate_buckling_resistance()
# The final answer is wrapped according to the format requirement.
# Example: <<<8177.9>>>
# I will use one decimal place as is common in engineering results.
print(f"<<<{final_answer:.1f}>>>")
