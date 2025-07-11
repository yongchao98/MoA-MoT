import math

def solve_moment_capacity():
    """
    Calculates the moment at collapse for the given reinforced concrete section.
    """
    # --- 1. Constants and Input Data ---
    fck = 30.0  # MPa
    fyk = 500.0 # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    Es = 200000.0 # MPa
    pi = math.pi

    # Geometry (in mm)
    d1 = 50.0   # Depth to top steel layer
    d2 = 260.0  # Depth to middle steel layer
    d3 = 350.0  # Depth to bottom steel layer
    h_trap = 300.0
    b_top = 100.0
    b_bottom = 400.0

    # Reinforcement
    D_bar = 20.0
    As_bar = pi * D_bar**2 / 4.0
    As1 = 2 * As_bar
    As2 = 2 * As_bar
    As3 = 3 * As_bar

    # --- 2. Material Design Properties ---
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    lambda_s = 0.8
    eta = 1.0
    stress_concrete = eta * fcd
    eps_cu = 0.0035
    eps_yd = fyd / Es

    # --- 3. Find Neutral Axis Depth 'x' by solving force equilibrium ---
    def calculate_net_force(x):
        """Calculates Net Force = Compression - Tension for a given x."""
        # Steel strains
        eps_s1 = eps_cu * (x - d1) / x if x > d1 else -eps_cu * (d1-x)/x
        eps_s2 = eps_cu * (d2 - x) / x
        eps_s3 = eps_cu * (d3 - x) / x

        # Steel stresses (positive for compression, negative for tension)
        fs1 = max(-fyd, min(fyd, eps_s1 * Es))
        fs2 = max(-fyd, min(fyd, eps_s2 * Es))
        fs3 = max(-fyd, min(fyd, eps_s3 * Es))

        # Concrete compression force
        s = lambda_s * x
        Fc_gross = 0
        if s > 0:
            if s <= h_trap:
                # Compression zone is a trapezoid. Width b(y) = 100 + y
                Ac_comp = 100 * s + s**2 / 2.0
                Fc_gross = stress_concrete * Ac_comp
            else:
                # This case is not expected here, but included for completeness
                Ac_trap = 100 * h_trap + h_trap**2 / 2.0
                Ac_rect = b_bottom * (s - h_trap)
                Fc_gross = stress_concrete * (Ac_trap + Ac_rect)

        # Steel forces (Positive = Compression, Negative = Tension)
        F_s1 = As1 * fs1
        F_s2 = As2 * fs2
        F_s3 = As3 * fs3
        
        # Total Compression Force
        # Account for displaced concrete by compression steel
        F_comp_total = 0
        if fs1 > 0: # If top steel is in compression
             F_comp_total += As1 * (fs1 - stress_concrete)
        else: # If top steel is in tension
             F_comp_total += F_s1

        if fs2 < 0: # Middle steel in tension
            F_tens_total = -F_s2
        else: # Middle steel in compression
            F_comp_total += As2 * (fs2 - stress_concrete)

        if fs3 < 0: # Bottom steel in tension
            F_tens_total = -F_s3
        else: # Bottom steel in compression
            F_comp_total += As3 * (fs3 - stress_concrete)

        
        # Net Force = Compression - Tension
        net_force = Fc_gross + F_comp_total - F_tens_total
        return net_force

    # Solve for x using bisection method
    x_low, x_high = d1, h_trap
    for _ in range(100):
        x_mid = (x_low + x_high) / 2
        if calculate_net_force(x_mid) > 0:
            x_high = x_mid
        else:
            x_low = x_mid
    x = (x_low + x_high) / 2
    
    # --- 4. Calculate Final Forces and Centroids with found x ---
    s = lambda_s * x
    
    # Concrete Force and Centroid
    Ac_comp = 100 * s + s**2 / 2
    Fc_gross = stress_concrete * Ac_comp
    b_s = b_top + s # Width at bottom of stress block
    zc = (s/3.0) * (b_top + 2*b_s) / (b_top + b_s)

    # Steel Strains and Stresses
    eps_s1 = eps_cu * (x - d1) / x
    eps_s2 = eps_cu * (d2 - x) / x
    eps_s3 = eps_cu * (d3 - x) / x
    fs1 = min(eps_s1 * Es, fyd)
    fs2 = min(eps_s2 * Es, fyd)
    fs3 = min(eps_s3 * Es, fyd)
    
    # Individual Forces (Magnitudes)
    Fsc1_net = As1 * (fs1 - stress_concrete)
    Fst2 = As2 * fs2
    Fst3 = As3 * fs3

    # Total Compression and Tension Forces
    Total_Compression = Fc_gross + Fsc1_net
    Total_Tension = Fst2 + Fst3
    
    # Centroids of total forces
    dc = (Fc_gross * zc + Fsc1_net * d1) / Total_Compression
    dt = (Fst2 * d2 + Fst3 * d3) / Total_Tension

    # --- 5. Calculate Ultimate Moment Mu ---
    lever_arm = dt - dc
    Mu_N_mm = Total_Tension * lever_arm
    Mu_kNm = Mu_N_mm / 1e6

    # --- Final Output ---
    print(f"CALCULATION STEPS:")
    print(f"1. Solved for neutral axis depth, x = {x:.2f} mm")
    print(f"2. Calculated total tensile force T = F_st_layer2 + F_st_layer3")
    print(f"   T = {Fst2/1000:.2f} kN + {Fst3/1000:.2f} kN = {Total_Tension/1000:.2f} kN")
    print(f"3. Calculated total compressive force C = F_concrete + F_sc_layer1_net")
    print(f"   C = {Fc_gross/1000:.2f} kN + {Fsc1_net/1000:.2f} kN = {Total_Compression/1000:.2f} kN (check: T ≈ C)")
    print(f"4. Calculated lever arm z = (centroid of tension) - (centroid of compression)")
    print(f"   z = {dt:.2f} mm - {dc:.2f} mm = {lever_arm:.2f} mm")
    print("\nFINAL EQUATION:")
    print(f"Moment at collapse, Mu = Total Tensile Force * Lever Arm")
    print(f"Mu = {Total_Tension:.2f} N * {lever_arm:.2f} mm")
    print(f"Mu = {Mu_N_mm:.2f} N.mm")
    print(f"Mu = {Mu_kNm:.2f} kNm")

solve_moment_capacity()