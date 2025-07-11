import math

def solve_collapse_moment():
    """
    Calculates the moment at collapse for a non-standard reinforced concrete section.
    """
    # 1. Define Material and Geometric Properties
    fck = 30.0  # MPa
    fyk = 500.0  # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    lambda_val = 0.8  # Factor for depth of rectangular stress block
    Es = 200000.0  # Modulus of Elasticity of steel in MPa
    epsilon_cu = 0.0035  # Ultimate concrete strain

    # Geometry (in mm)
    b_top = 100.0
    b_bottom = 400.0
    h_trap = 300.0
    h_rect = 100.0
    h_total = h_trap + h_rect

    # Reinforcement
    d_s1 = 50.0   # Depth of top steel (2H20)
    d_s2 = 50.0 + 210.0  # Depth of middle steel (2H20)
    d_s3 = 50.0 + 210.0 + 90.0 # Depth of bottom steel (3H20)
    
    dia_bar = 20.0
    area_one_bar = math.pi * (dia_bar / 2)**2
    
    As1 = 2 * area_one_bar # Compression steel
    As2 = 2 * area_one_bar # Tension steel
    As3 = 3 * area_one_bar # Tension steel

    # 2. Calculate Design Strengths
    f_cd = alpha_cc * fck / gamma_c
    f_yd = fyk / gamma_s
    epsilon_yd = f_yd / Es

    print("--- Design Parameters ---")
    print(f"f_cd (Design concrete strength): {f_cd:.2f} MPa")
    print(f"f_yd (Design steel yield strength): {f_yd:.2f} MPa")
    print(f"Steel yield strain: {epsilon_yd:.5f}")
    print(f"Area of top steel (As1): {As1:.2f} mm^2")
    print(f"Area of middle steel (As2): {As2:.2f} mm^2")
    print(f"Area of bottom steel (As3): {As3:.2f} mm^2")
    print("-" * 25)

    # 3. Locate the Neutral Axis (x) using Bisection Method
    low_x = 0.0
    high_x = h_total
    x = 0
    for _ in range(100): # 100 iterations for high precision
        x = (low_x + high_x) / 2
        s = lambda_val * x # Depth of stress block

        # --- Calculate Forces for current x ---
        
        # Concrete Force (Fc) and its centroid (yc)
        # Assuming s is within the trapezoidal part (s <= 300 mm)
        # This is a reasonable assumption that will be verified by the final x value
        if s <= h_trap:
            b_at_s = b_top + s # Width at depth s
            Ac = (b_top + b_at_s) / 2 * s
            Fc = f_cd * Ac
            # Centroid of the trapezoidal compression block from the top fiber
            yc = (s / 3) * (b_top + 2 * b_at_s) / (b_top + b_at_s)
        else: # Full trapezoid + a rectangle
             Ac_trap = (b_top + b_bottom) / 2 * h_trap
             yc_trap = (h_trap / 3) * (b_top + 2 * b_bottom) / (b_top + b_bottom)
             
             s_rect = s - h_trap
             Ac_rect = b_bottom * s_rect
             yc_rect = h_trap + s_rect / 2
             
             Ac = Ac_trap + Ac_rect
             Fc = f_cd * Ac
             yc = (Ac_trap * yc_trap + Ac_rect * yc_rect) / Ac

        # Steel Forces
        # Top steel (compression)
        epsilon_s1 = epsilon_cu * (x - d_s1) / x if x > 0 else 0
        sigma_s1 = max(-f_yd, min(f_yd, Es * epsilon_s1))
        # Deduct force of displaced concrete
        Fsc = As1 * (sigma_s1 - f_cd) if s > d_s1 else As1 * sigma_s1
        
        # Middle steel (tension)
        epsilon_s2 = epsilon_cu * (d_s2 - x) / x if x > 0 else 0
        sigma_s2 = max(-f_yd, min(f_yd, Es * epsilon_s2))
        Fst2 = As2 * sigma_s2
        
        # Bottom steel (tension)
        epsilon_s3 = epsilon_cu * (d_s3 - x) / x if x > 0 else 0
        sigma_s3 = max(-f_yd, min(f_yd, Es * epsilon_s3))
        Fst3 = As3 * sigma_s3
        
        # Check force equilibrium
        total_force = Fc + Fsc - Fst2 - Fst3
        
        if total_force > 0: # Too much compression, neutral axis is too deep
            high_x = x
        else:
            low_x = x

    print("--- Equilibrium Calculation ---")
    print(f"Final neutral axis depth (x): {x:.2f} mm")
    s = lambda_val * x
    print(f"Final stress block depth (s = 0.8x): {s:.2f} mm")

    # 4. Calculate Final Forces and Collapse Moment (Mu)
    # Recalculate all forces and centroids with the final 'x'
    # Concrete Force (Fc) and its centroid (yc)
    b_at_s = b_top + s
    Ac = (b_top + b_at_s) / 2 * s
    Fc = f_cd * Ac
    yc = (s / 3) * (b_top + 2 * b_at_s) / (b_top + b_at_s)

    # Steel Forces
    epsilon_s1 = epsilon_cu * (x - d_s1) / x
    sigma_s1 = max(-f_yd, min(f_yd, Es * epsilon_s1))
    Fsc = As1 * (sigma_s1 - f_cd)
    
    epsilon_s2 = epsilon_cu * (d_s2 - x) / x
    sigma_s2 = max(-f_yd, min(f_yd, Es * epsilon_s2))
    Fst2 = As2 * sigma_s2

    epsilon_s3 = epsilon_cu * (d_s3 - x) / x
    sigma_s3 = max(-f_yd, min(f_yd, Es * epsilon_s3))
    Fst3 = As3 * sigma_s3
    
    print("\n--- Final Forces at Collapse ---")
    print(f"Concrete compressive force (Fc): {Fc/1000:.2f} kN")
    print(f"Compression steel force (Fsc): {Fsc/1000:.2f} kN")
    print(f"Middle tension steel force (Fst2): {Fst2/1000:.2f} kN (Stress: {sigma_s2:.2f} MPa, Strain: {epsilon_s2:.5f})")
    print(f"Bottom tension steel force (Fst3): {Fst3/1000:.2f} kN (Stress: {sigma_s3:.2f} MPa, Strain: {epsilon_s3:.5f})")
    Total_Comp = (Fc + Fsc) / 1000
    Total_Tens = (Fst2 + Fst3) / 1000
    print(f"Total Compression = {Total_Comp:.2f} kN, Total Tension = {Total_Tens:.2f} kN")
    
    # Calculate moment about the top fiber (in N.mm)
    Mc = Fc * yc
    Msc = Fsc * d_s1
    Mst2 = Fst2 * d_s2
    Mst3 = Fst3 * d_s3
    
    Mu_Nmm = Mst2 + Mst3 - Msc - Mc
    Mu_kNm = Mu_Nmm / 1e6

    print("\n--- Collapse Moment Calculation (about top fiber) ---")
    print(f"M_u = M_st_mid + M_st_bot - M_sc - M_c")
    print(f"M_u = ({Fst2/1000:.2f} kN * {d_s2} mm) + ({Fst3/1000:.2f} kN * {d_s3} mm) - ({Fsc/1000:.2f} kN * {d_s1} mm) - ({Fc/1000:.2f} kN * {yc:.2f} mm)")
    print(f"M_u = {Mst2/1e6:.2f} + {Mst3/1e6:.2f} - {Msc/1e6:.2f} - {Mc/1e6:.2f} kNm")
    print(f"M_u = {Mu_kNm:.2f} kNm")
    
    # Final answer
    print("\nFinal Answer:")
    print(f"The moment at collapse is {Mu_kNm:.1f} kNm.")
    print(f"<<<{Mu_kNm:.1f}>>>")

solve_collapse_moment()