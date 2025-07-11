import math

def solve_buckling_resistance():
    """
    Calculates the buckling resistance of a concrete encased composite column
    based on Eurocode 4 (EN 1994-1-1).
    """
    # Step 1: Define Material and Section Properties
    print("--- Step 1: Defining Input Properties ---")
    
    # Column geometry
    L = 4000  # Storey height in mm
    H = 400   # Column height in mm
    B = 400   # Column width in mm
    cover = 30 # Concrete cover in mm

    # Structural Steel (UC 254x254x132, S355)
    f_y = 345      # Yield strength in N/mm^2
    A_a = 16800    # Area in mm^2 (168 cm^2)
    I_az = 7530e4  # Second moment of area about z-axis in mm^4 (7530 cm^4)
    E_a = 210000   # Modulus of elasticity for steel in N/mm^2

    # Concrete (C40/50)
    f_ck = 40      # Characteristic cylinder strength in N/mm^2

    # Reinforcing Steel (4 x T16, S500)
    rebar_dia = 16 # Rebar diameter in mm
    num_rebars = 4
    f_sk = 500     # Characteristic yield strength of rebar in N/mm^2
    E_s = 200000   # Modulus of elasticity for rebar in N/mm^2

    # Partial Safety Factors (EN 1994-1-1)
    gamma_a = 1.1  # For structural steel in compression
    gamma_c = 1.5  # For concrete
    gamma_s = 1.15 # For reinforcing steel

    print(f"Column Length (L): {L} mm")
    print(f"Steel Yield Strength (f_y): {f_y} N/mm^2")
    print(f"Concrete Strength (f_ck): {f_ck} N/mm^2")
    print(f"Rebar Strength (f_sk): {f_sk} N/mm^2\n")

    # Step 2: Calculate Additional Sectional Properties
    print("--- Step 2: Calculating Sectional Properties ---")
    
    # Reinforcing steel area
    A_s_one_bar = math.pi * (rebar_dia / 2)**2
    A_s = num_rebars * A_s_one_bar
    print(f"Total Rebar Area (A_s): {A_s:.2f} mm^2")

    # Concrete area (net)
    A_c_gross = H * B
    A_c = A_c_gross - A_a - A_s
    print(f"Net Concrete Area (A_c): {A_c:.2f} mm^2")

    # Second moment of area for rebars about z-axis
    # Rebars are in the corners. Distance from z-axis (horizontal centerline) to rebar centroid.
    d_rebar_y = H / 2 - cover - rebar_dia / 2
    I_sz = num_rebars * A_s_one_bar * d_rebar_y**2
    # Note: We are checking buckling about the z-axis (weaker axis of I-beam). 
    # The distance from the z-axis to the rebars is d_rebar_y.
    print(f"Distance of rebar from z-axis: {d_rebar_y} mm")
    print(f"Rebar Moment of Inertia (I_sz): {I_sz:.2e} mm^4")

    # Second moment of area for concrete about z-axis
    I_c_gross = B * H**3 / 12
    I_cz = I_c_gross - I_az 
    print(f"Concrete Moment of Inertia (I_cz): {I_cz:.2e} mm^4\n")

    # Step 3: Calculate Effective Flexural Stiffness ((EI)_eff)
    print("--- Step 3: Calculating Effective Flexural Stiffness (EI)_eff ---")
    # Secant modulus of elasticity of concrete
    E_cm = 22 * ((f_ck + 8) / 10)**(1/3) * 1000 # in N/mm^2
    print(f"Concrete Modulus of Elasticity (E_cm): {E_cm:.0f} N/mm^2")
    
    # K_e = 1.0 because the problem states to ignore creep and shrinkage
    K_e = 1.0
    EI_eff_z = E_a * I_az + K_e * E_cm * I_cz + E_s * I_sz
    print(f"Effective Flexural Stiffness ((EI)_eff,z): {EI_eff_z:.3e} Nmm^2\n")

    # Step 4: Calculate Design Plastic Resistance (N_pl,Rd)
    print("--- Step 4: Calculating Design Plastic Resistance (N_pl,Rd) ---")
    f_yd = f_y / gamma_a
    f_cd = f_ck / gamma_c
    f_sd = f_sk / gamma_s
    alpha_c = 0.85 # for encased sections
    
    N_pl_Rd = A_a * f_yd + alpha_c * A_c * f_cd + A_s * f_sd
    print(f"N_pl,Rd = A_a*f_yd + α_c*A_c*f_cd + A_s*f_sd")
    print(f"N_pl,Rd = {A_a:.0f}*{f_yd:.2f} + {alpha_c:.2f}*{A_c:.0f}*{f_cd:.2f} + {A_s:.0f}*{f_sd:.2f}")
    print(f"N_pl,Rd = {A_a*f_yd/1000:.1f} kN + {alpha_c*A_c*f_cd/1000:.1f} kN + {A_s*f_sd/1000:.1f} kN")
    print(f"Design Plastic Resistance (N_pl,Rd): {N_pl_Rd / 1000:.1f} kN\n")

    # Step 5: Calculate Critical Buckling Load (N_cr)
    print("--- Step 5: Calculating Critical Buckling Load (N_cr) ---")
    # Assuming pinned-pinned ends, effective length L_cr = 1.0 * L
    L_cr = 1.0 * L
    N_cr_z = (math.pi**2 * EI_eff_z) / L_cr**2
    print(f"Effective Length (L_cr): {L_cr:.0f} mm")
    print(f"Critical Buckling Load (N_cr): {N_cr_z / 1000:.1f} kN\n")

    # Step 6: Determine Relative Slenderness (lambda_bar)
    print("--- Step 6: Determining Relative Slenderness (λ_bar) ---")
    # Characteristic plastic resistance (using characteristic strengths)
    N_pl_Rk = A_a * f_y + alpha_c * A_c * f_ck + A_s * f_sk
    print(f"Characteristic Plastic Resistance (N_pl,Rk): {N_pl_Rk / 1000:.1f} kN")
    
    lambda_bar = math.sqrt(N_pl_Rk / N_cr_z)
    print(f"λ_bar = sqrt(N_pl,Rk / N_cr) = sqrt({N_pl_Rk/1000:.1f} / {N_cr_z/1000:.1f})")
    print(f"Relative Slenderness (λ_bar): {lambda_bar:.3f}\n")

    # Step 7: Find Buckling Reduction Factor (chi)
    print("--- Step 7: Calculating Buckling Reduction Factor (χ) ---")
    # For encased I-section, buckling about z-z axis -> Buckling Curve 'c'
    alpha_imp = 0.49 # Imperfection factor for curve 'c'
    
    phi = 0.5 * (1 + alpha_imp * (lambda_bar - 0.2) + lambda_bar**2)
    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    
    print(f"Using Buckling Curve 'c' (α = {alpha_imp})")
    print(f"Φ = 0.5 * [1 + α*(λ_bar-0.2) + λ_bar^2] = {phi:.3f}")
    print(f"χ = 1 / (Φ + sqrt(Φ^2 - λ_bar^2)) = 1 / ({phi:.3f} + sqrt({phi:.3f}^2 - {lambda_bar:.3f}^2))")
    print(f"Buckling Reduction Factor (χ): {chi:.3f}\n")

    # Step 8: Calculate Buckling Resistance (N_b,Rd)
    print("--- Step 8: Final Buckling Resistance (N_b,Rd) ---")
    N_b_Rd = chi * N_pl_Rd
    
    print(f"N_b,Rd = χ * N_pl,Rd")
    print(f"N_b,Rd = {chi:.3f} * {N_pl_Rd / 1000:.1f} kN")
    print(f"The final buckling resistance is: {N_b_Rd / 1000:.1f} kN")
    
    return N_b_Rd / 1000

# Execute the function and print the final answer in the required format
final_answer = solve_buckling_resistance()
print(f"\n<<< {final_answer:.1f} >>>")