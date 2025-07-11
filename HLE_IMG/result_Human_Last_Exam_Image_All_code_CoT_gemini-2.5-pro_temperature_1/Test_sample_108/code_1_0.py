import math

def solve_moment_at_collapse():
    """
    Calculates the moment at collapse for the given non-standard concrete section.
    """

    # 1. Define Material Properties and Parameters
    fck = 30.0  # MPa
    fyk = 500.0  # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    Es = 200000.0  # MPa
    epsilon_cu3 = 0.0035
    lambda_param = 0.8

    # Design strengths
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    epsilon_yd = fyd / Es

    print("--- Material Properties ---")
    print(f"Design compressive strength of concrete, fcd = {alpha_cc:.2f} * {fck} / {gamma_c} = {fcd:.2f} MPa")
    print(f"Design yield strength of steel, fyd = {fyk} / {gamma_s} = {fyd:.2f} MPa")
    print("-" * 20)

    # 2. Define Geometry and Reinforcement
    # Section geometry (all units in mm)
    # Origin (y_top = 0) is at the top fiber
    # Width b(y_top) = y_top + 100 for 0 <= y_top <= 300
    # Width b(y_top) = 400 for 300 < y_top <= 400

    # Reinforcement
    D_bar = 20.0  # mm
    As_bar = math.pi * (D_bar ** 2) / 4

    # Layer 1 (Compression)
    d1 = 50.0
    As1 = 2 * As_bar

    # Layer 2 (Tension)
    d2 = 50.0 + 210.0
    As2 = 2 * As_bar

    # Layer 3 (Tension)
    d3 = d2 + 90.0
    As3 = 3 * As_bar
    
    print("--- Reinforcement Details ---")
    print(f"Area of one H20 bar: {As_bar:.2f} mm^2")
    print(f"Layer 1 (Compression): 2H20 at d1 = {d1} mm, Area As1 = {As1:.2f} mm^2")
    print(f"Layer 2 (Tension): 2H20 at d2 = {d2} mm, Area As2 = {As2:.2f} mm^2")
    print(f"Layer 3 (Tension): 3H20 at d3 = {d3} mm, Area As3 = {As3:.2f} mm^2")
    print("-" * 20)


    # 3. Find Neutral Axis Depth (x) using bisection method
    
    def get_net_force(x):
        """Calculates Compression Force - Tension Force for a given neutral axis depth x."""
        # --- Compression Forces ---
        
        # Concrete compression force (Fc)
        s = lambda_param * x  # Depth of rectangular stress block
        if s <= 300:
            # Compression zone is within the top trapezoid
            Ac = s**2 / 2 + 100 * s
        else:
            # Compression zone extends into the rectangular part
            Ac_trap = 300**2 / 2 + 100 * 300 # Area of the top 300mm trapezoid part
            Ac_rect = 400 * (s - 300)
            Ac = Ac_trap + Ac_rect
        Fc = fcd * Ac

        # Steel compression force (Fsc)
        # Strain in compression steel
        epsilon_sc = epsilon_cu3 * (x - d1) / x if x > d1 else 0
        # Stress in compression steel (accounting for yield)
        sigma_sc = min(epsilon_sc * Es, fyd)
        # Force in compression steel (ignoring displaced concrete as is common practice)
        Fsc = As1 * sigma_sc

        # Total Compression Force
        C_total = Fc + Fsc

        # --- Tension Forces ---
        
        # Steel tension force in layer 2 (Fst2)
        if x < d2:
            epsilon_s2 = epsilon_cu3 * (d2 - x) / x
            sigma_s2 = min(epsilon_s2 * Es, fyd)
            Fst2 = As2 * sigma_s2
        else: # Layer 2 is in compression
            Fst2 = 0 
        
        # Steel tension force in layer 3 (Fst3)
        if x < d3:
            epsilon_s3 = epsilon_cu3 * (d3 - x) / x
            sigma_s3 = min(epsilon_s3 * Es, fyd)
            Fst3 = As3 * sigma_s3
        else: # Layer 3 is in compression
            Fst3 = 0

        # Total Tension Force
        T_total = Fst2 + Fst3
        
        return C_total - T_total

    # Bisection search
    low = d1
    high = d3
    for _ in range(100): # 100 iterations for high precision
        x_guess = (low + high) / 2
        net_force = get_net_force(x_guess)
        if net_force > 0: # Too much compression, x is too large
            high = x_guess
        else: # Too little compression, x is too small
            low = x_guess
    
    x = (low + high) / 2
    print("--- Equilibrium Calculation ---")
    print(f"Calculated neutral axis depth, x = {x:.2f} mm")
    print("-" * 20)
    
    # 4. Calculate Moment at Collapse (M_Rd)
    
    # Recalculate forces with final x
    s = lambda_param * x
    if s <= 300:
        Ac = s**2 / 2 + 100 * s
        moment_of_area = s**3 / 3 + 50 * s**2 # Moment of area about top fiber
        z_c = moment_of_area / Ac if Ac > 0 else 0
    else:
        Ac_trap = 300**2 / 2 + 100 * 300
        moment_trap = 300**3 / 3 + 50 * 300**2
        
        Ac_rect = 400 * (s - 300)
        # Moment of rect area about top fiber = integral from 300 to s of y*400 dy
        moment_rect = 200 * (s**2 - 300**2) 
        
        Ac = Ac_trap + Ac_rect
        moment_of_area = moment_trap + moment_rect
        z_c = moment_of_area / Ac if Ac > 0 else 0

    Fc = fcd * Ac
    
    epsilon_sc = epsilon_cu3 * (x - d1) / x if x > d1 else 0
    sigma_sc = min(epsilon_sc * Es, fyd)
    Fsc = As1 * sigma_sc
    
    epsilon_s2 = epsilon_cu3 * (d2 - x) / x
    sigma_s2 = min(epsilon_s2 * Es, fyd)
    Fst2 = As2 * sigma_s2
    
    epsilon_s3 = epsilon_cu3 * (d3 - x) / x
    sigma_s3 = min(epsilon_s3 * Es, fyd)
    Fst3 = As3 * sigma_s3

    # Moment calculation about the top fiber
    # M_Rd = Sum of (Force * lever arm)
    M_Rd_Nmm = (Fst2 * d2) + (Fst3 * d3) - (Fsc * d1) - (Fc * z_c)
    M_Rd_kNm = M_Rd_Nmm / 1e6
    
    print("--- Final Forces and Moment Calculation ---")
    print(f"Concrete compressive force, Fc = {Fc/1000:.2f} kN")
    print(f"Centroid of concrete compression block, z_c = {z_c:.2f} mm")
    print(f"Compression steel force, Fsc = {Fsc/1000:.2f} kN")
    print(f"Tension steel force (Layer 2), Fst2 = {Fst2/1000:.2f} kN")
    print(f"Tension steel force (Layer 3), Fst3 = {Fst3/1000:.2f} kN")
    
    total_compression = (Fc + Fsc) / 1000
    total_tension = (Fst2 + Fst3) / 1000
    print(f"Total Compression = {total_compression:.2f} kN, Total Tension = {total_tension:.2f} kN (check)")
    print("-" * 20)

    print("Moment at collapse, M_Rd = Fst2*d2 + Fst3*d3 - Fsc*d1 - Fc*z_c")
    print(f"M_Rd = ({Fst2/1000:.2f} kN * {d2:.1f} mm) + ({Fst3/1000:.2f} kN * {d3:.1f} mm) - ({Fsc/1000:.2f} kN * {d1:.1f} mm) - ({Fc/1000:.2f} kN * {z_c:.2f} mm)")
    
    M_tens_kNm = ((Fst2 * d2) + (Fst3 * d3)) / 1e6
    M_comp_kNm = ((Fsc * d1) + (Fc * z_c)) / 1e6
    print(f"M_Rd = {M_tens_kNm:.2f} kNm - {M_comp_kNm:.2f} kNm")

    print(f"\nFinal Moment at Collapse, M_Rd = {M_Rd_kNm:.2f} kNm")

    return M_Rd_kNm

final_answer = solve_moment_at_collapse()
print(f'<<<{final_answer:.2f}>>>')