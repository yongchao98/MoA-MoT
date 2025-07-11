import math

def solve_moment_at_collapse():
    """
    Calculates the moment at collapse for the given non-standard RC section.
    """
    # --- 1. Define Constants and Material Properties ---
    print("--- 1. Material Properties and Design Strengths ---")
    fck = 30.0  # MPa
    fyk = 500.0 # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    eta = 1.0
    lambda_ = 0.8
    E_s = 200000.0 # MPa, Modulus of Elasticity of steel

    # Design strengths
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    eps_cu3 = 0.0035  # Ultimate compressive strain in concrete
    eps_yd = fyd / E_s    # Yield strain of steel

    print(f"Design compressive strength of concrete, f_cd = {fcd:.2f} MPa")
    print(f"Design yield strength of steel, f_yd = {fyd:.2f} MPa")
    print(f"Yield strain of steel, eps_yd = {eps_yd:.5f}")
    print("-" * 50)

    # --- 2. Section Geometry and Steel Areas ---
    print("--- 2. Section Geometry and Steel Areas ---")
    # Geometry (all in mm)
    d_top = 50.0   # Depth of top (compression) steel
    d_mid = 260.0  # Depth of middle (tension) steel
    d_bot = 350.0  # Depth of bottom (tension) steel
    h_trap = 300.0 # Height of the upper trapezoidal part
    w_top = 100.0  # Width at the top of the section

    # Steel Areas (H20 bar diameter = 20mm)
    d_bar = 20.0
    A_s_one_bar = math.pi * (d_bar / 2)**2
    A_s_top = 2 * A_s_one_bar
    A_s_mid = 2 * A_s_one_bar
    A_s_bot = 3 * A_s_one_bar

    print(f"Area of top compression steel (2H20), A_s' = {A_s_top:.2f} mm^2")
    print(f"Area of middle tension steel (2H20), A_s1 = {A_s_mid:.2f} mm^2")
    print(f"Area of bottom tension steel (3H20), A_s2 = {A_s_bot:.2f} mm^2")
    print("-" * 50)
    
    # --- 3. Iterative Search for Neutral Axis (x) ---
    print("--- 3. Finding the Neutral Axis Depth (x) ---")
    print("Solving C = T by iterating for x, where C = F_cc + F_sc and T = F_st,mid + F_st,bot")

    def calculate_forces(x):
        # Depth of equivalent rectangular stress block
        s = lambda_ * x
        
        # a) Concrete compression force (F_cc)
        # The section width b(y) = 100 + y for y <= 300 mm
        A_cc = 100 * s + 0.5 * s**2
        F_cc = eta * fcd * A_cc

        # b) Compression steel force (F_sc)
        F_sc = 0
        if x > d_top:
            eps_sc = eps_cu3 * (x - d_top) / x
            sigma_sc = min(eps_sc * E_s, fyd)
            # Force in steel minus force in displaced concrete
            F_sc = A_s_top * (sigma_sc - eta * fcd)
            F_sc = max(0, F_sc) # Force cannot be tensile in the compression force sum
        
        # c) Middle tension steel force (F_st_mid)
        F_st_mid = 0
        if x < d_mid:
            eps_st_mid = eps_cu3 * (d_mid - x) / x
            sigma_st_mid = min(eps_st_mid * E_s, fyd)
            F_st_mid = A_s_mid * sigma_st_mid

        # d) Bottom tension steel force (F_st_bot)
        F_st_bot = 0
        if x < d_bot:
            eps_st_bot = eps_cu3 * (d_bot - x) / x
            sigma_st_bot = min(eps_st_bot * E_s, fyd)
            F_st_bot = A_s_bot * sigma_st_bot
            
        Total_C = F_cc + F_sc
        Total_T = F_st_mid + F_st_bot
        
        return Total_C, Total_T

    # Bisection method to find x
    x_low, x_high = 1.0, h_trap / lambda_ # Initial search range
    for _ in range(100): # 100 iterations for high precision
        x_guess = (x_low + x_high) / 2
        C, T = calculate_forces(x_guess)
        if C > T:
            x_high = x_guess
        else:
            x_low = x_guess
    x = (x_low + x_high) / 2
    
    print(f"Equilibrium found at neutral axis depth, x = {x:.2f} mm")
    print("-" * 50)

    # --- 4. Final Force Calculation ---
    print("--- 4. Calculating Final Internal Forces ---")
    s = lambda_ * x
    A_cc = 100 * s + 0.5 * s**2
    F_cc = eta * fcd * A_cc

    eps_sc = eps_cu3 * (x - d_top) / x
    sigma_sc = min(eps_sc * E_s, fyd)
    F_sc = max(0, A_s_top * (sigma_sc - eta * fcd))

    eps_st_mid = eps_cu3 * (d_mid - x) / x
    sigma_st_mid = min(eps_st_mid * E_s, fyd)
    F_st_mid = A_s_mid * sigma_st_mid

    eps_st_bot = eps_cu3 * (d_bot - x) / x
    sigma_st_bot = min(eps_st_bot * E_s, fyd)
    F_st_bot = A_s_bot * sigma_st_bot
    
    Total_C = F_cc + F_sc
    Total_T = F_st_mid + F_st_bot

    print(f"Force in concrete compression block, F_cc = {F_cc/1000:.2f} kN")
    print(f"Force in compression steel, F_sc = {F_sc/1000:.2f} kN")
    print(f"Force in middle tension steel, F_st,mid = {F_st_mid/1000:.2f} kN")
    print(f"Force in bottom tension steel, F_st,bot = {F_st_bot/1000:.2f} kN")
    print(f"Total Compression = {Total_C/1000:.2f} kN | Total Tension = {Total_T/1000:.2f} kN")
    print("-" * 50)

    # --- 5. Calculate Moment of Resistance (M_Rd) ---
    print("--- 5. Calculating Moment at Collapse (M_Rd) ---")
    # Centroid of the trapezoidal concrete compression area from the top fiber
    # Formula for centroid of a trapezoid from its top base 'a'
    # a = 100, b = 100+s, h = s
    z_cc = (s / 3.0) * (2 * (100 + s) + 100) / (100 + (100 + s))
    
    print(f"Lever arm for F_cc (z_cc) = {z_cc:.2f} mm")
    print(f"Lever arm for F_sc (z_sc) = {d_top:.2f} mm")
    print(f"Lever arm for F_st,mid (z_s1) = {d_mid:.2f} mm")
    print(f"Lever arm for F_st,bot (z_s2) = {d_bot:.2f} mm")
    print("\nCalculating moment about the top fiber:")
    print("M_Rd = (F_st,mid * z_s1) + (F_st,bot * z_s2) - (F_cc * z_cc) - (F_sc * z_sc)")

    M_Rd = (F_st_mid * d_mid) + (F_st_bot * d_bot) - (F_cc * z_cc) - (F_sc * d_top)
    M_Rd_kNm = M_Rd / 1e6
    
    # Print the equation with final numbers
    print(f"M_Rd = ({F_st_mid/1000:.2f} kN * {d_mid/1000:.3f} m) + ({F_st_bot/1000:.2f} kN * {d_bot/1000:.3f} m) - ({F_cc/1000:.2f} kN * {z_cc/1000:.3f} m) - ({F_sc/1000:.2f} kN * {d_top/1000:.3f} m)")

    print("\n" + "="*50)
    print(f"The moment at collapse is: {M_Rd_kNm:.2f} kNm")
    print("="*50)
    
    return M_Rd_kNm

# Run the solver and print the final result
final_moment = solve_moment_at_collapse()
print(f"<<<{final_moment:.2f}>>>")
