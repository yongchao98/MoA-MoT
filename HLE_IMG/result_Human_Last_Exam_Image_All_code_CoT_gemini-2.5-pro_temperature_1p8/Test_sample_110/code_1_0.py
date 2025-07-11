import math

def calculate_buckling_resistance():
    """
    Calculates the buckling resistance of a concrete encased composite column.
    All units are in N and mm unless specified otherwise.
    """
    # 1. Define Material and Geometric Properties from the problem statement
    print("--- 1. Material and Geometric Properties ---")
    
    # Steel Section (UC 254x254x132 S355)
    f_y = 345.0  # N/mm^2, characteristic yield strength
    A_a = 168.0 * 100  # mm^2 (168 cm^2)
    I_a_y = 22500.0 * 1e4 # mm^4 (strong axis)
    I_a_z = 7530.0 * 1e4  # mm^4 (weak axis)
    E_a = 210000.0  # N/mm^2, Modulus of Elasticity of structural steel
    
    print(f"Structural Steel (S355): fy = {f_y} MPa, A_a = {A_a} mm^2, I_y = {I_a_y:.2e} mm^4, I_z = {I_a_z:.2e} mm^4")

    # Concrete (C40/50)
    f_ck = 40.0  # N/mm^2, characteristic cylinder strength
    h_c_gross = 400.0  # mm
    b_c_gross = 400.0  # mm
    
    print(f"Concrete (C40/50): fck = {f_ck} MPa, Gross Section = {h_c_gross}x{b_c_gross} mm")
    
    # Reinforcing Steel (4 x T16 S500)
    d_rebar = 16.0  # mm
    num_rebars = 4
    f_sk = 500.0  # N/mm^2, characteristic yield strength
    E_s = 200000.0 # N/mm^2, Modulus of Elasticity of reinforcing steel
    cover = 30.0   # mm
    
    print(f"Reinforcement (S500): {num_rebars} no. T{int(d_rebar)} bars, fsk = {f_sk} MPa")
    
    # Column parameters
    L = 4.0 * 1000 # mm, storey height
    L_cr = L       # mm, effective length (assuming pinned-pinned)
    print(f"Column Length: L = {L} mm, Effective Length L_cr = {L_cr} mm")
    
    # 2. Calculate Sectional Properties
    print("\n--- 2. Sectional Property Calculations ---")
    
    # Net Concrete Area
    A_c = (h_c_gross * b_c_gross) - A_a
    print(f"Net Concrete Area (A_c): {A_c:.2f} mm^2")

    # Second Moment of Area of Concrete (net)
    I_c_gross = (b_c_gross * h_c_gross**3) / 12
    I_c_y = I_c_gross - I_a_y
    I_c_z = I_c_gross - I_a_z
    print(f"Net Concrete Second Moment of Area: I_c_y = {I_c_y:.2e} mm^4, I_c_z = {I_c_z:.2e} mm^4")
    
    # Reinforcement Area and Second Moment of Area
    A_s1 = math.pi * (d_rebar**2) / 4
    A_s = num_rebars * A_s1
    dist_rebar_center = h_c_gross/2 - cover - d_rebar/2 # Distance from section centroid to rebar centroid
    I_s = num_rebars * A_s1 * dist_rebar_center**2  # Same for y and z axes due to symmetry
    print(f"Reinforcement Area (A_s): {A_s:.2f} mm^2")
    print(f"Reinforcement Second Moment of Area (I_s): {I_s:.2e} mm^4 (for both axes)")

    # 3. Calculate Effective Flexural Stiffness (EI_eff)
    print("\n--- 3. Effective Flexural Stiffness (EI)_eff ---")

    # Secant Modulus of Concrete
    E_cm = 22 * ((f_ck + 8) / 10)**0.3
    print(f"Concrete Secant Modulus (E_cm): {E_cm:.2f} MPa")
    
    # Since problem states to ignore creep, Ke = 1.0
    K_e = 1.0
    print(f"Stiffness Factor for Concrete (K_e): {K_e} (assuming 'ignore creep' means using short-term properties)")

    EI_eff_y = E_a * I_a_y + E_s * I_s + K_e * E_cm * I_c_y
    EI_eff_z = E_a * I_a_z + E_s * I_s + K_e * E_cm * I_c_z
    
    print(f"(EI)_eff,y = {EI_eff_y:.3e} Nmm^2")
    print(f"(EI)_eff,z = {EI_eff_z:.3e} Nmm^2")

    # The governing stiffness is the smaller value
    EI_eff_gov = min(EI_eff_y, EI_eff_z)
    print(f"Governing stiffness is (EI)_eff_gov = {EI_eff_gov:.3e} Nmm^2 (buckling about z-z axis)")
    
    # 4. Determine Critical Buckling Load (N_cr)
    print("\n--- 4. Critical Buckling Load N_cr ---")
    N_cr = (math.pi**2 * EI_eff_gov) / (L_cr**2)
    print(f"N_cr = (pi^2 * {EI_eff_gov:.3e}) / {L_cr}^2 = {N_cr/1000:.2f} kN")
    
    # 5. Calculate Plastic Resistance (N_pl)
    print("\n--- 5. Plastic Resistance N_pl ---")
    
    # Characteristic Resistance (N_pl,Rk)
    N_pl_Rk = A_a * f_y + A_c * f_ck + A_s * f_sk
    print(f"Characteristic Resistance (N_pl,Rk) = {N_pl_Rk/1000:.2f} kN")

    # Design Resistance (N_pl,Rd)
    gamma_a = 1.0
    gamma_c = 1.5
    gamma_s = 1.15
    f_yd = f_y / gamma_a
    f_cd = f_ck / gamma_c  # Assuming alpha_cc = 1.0 as per EC2
    f_sd = f_sk / gamma_s
    
    # Per EC4 for encased sections, alpha_c = 1.0 for resistance calculation
    alpha_c = 1.0
    N_pl_Rd = A_a * f_yd + alpha_c * A_c * f_cd + A_s * f_sd
    print(f"Design Resistance (N_pl,Rd) = {N_pl_Rd/1000:.2f} kN")

    # 6. Calculate Slenderness and Buckling Factors
    print("\n--- 6. Slenderness and Buckling Factors ---")
    
    # Non-dimensional Slenderness (lambda_bar)
    lambda_bar = math.sqrt(N_pl_Rk / N_cr)
    print(f"Non-dimensional Slenderness (lambda_bar) = sqrt({N_pl_Rk/1000:.2f} / {N_cr/1000:.2f}) = {lambda_bar:.4f}")
    
    # Buckling Curve (z-z axis of encased I-section -> curve 'b')
    alpha_buckling = 0.34 # Imperfection factor for buckling curve 'b'
    print(f"Buckling curve 'b' is used (alpha = {alpha_buckling})")

    # Reduction factor (chi)
    phi = 0.5 * (1 + alpha_buckling * (lambda_bar - 0.2) + lambda_bar**2)
    chi = 1 / (phi + math.sqrt(phi**2 - lambda_bar**2))
    if chi > 1.0:
      chi = 1.0
    print(f"phi = {phi:.4f}")
    print(f"Buckling Reduction Factor (chi) = {chi:.4f}")

    # 7. Calculate Final Buckling Resistance (N_b,Rd)
    print("\n--- 7. Final Buckling Resistance N_b,Rd ---")
    
    N_b_Rd = chi * N_pl_Rd
    
    print(f"The design buckling resistance is given by N_b,Rd = chi * N_pl,Rd")
    print(f"N_b,Rd = {chi:.4f} * {N_pl_Rd/1000:.2f} kN")
    print(f"N_b,Rd = {N_b_Rd/1000:.1f} kN")
    
    return N_b_Rd / 1000

final_answer = calculate_buckling_resistance()
print(f"\nThe final buckling resistance of the column is {final_answer:.1f} kN.")
print(f"\n<<<{final_answer:.1f}>>>")