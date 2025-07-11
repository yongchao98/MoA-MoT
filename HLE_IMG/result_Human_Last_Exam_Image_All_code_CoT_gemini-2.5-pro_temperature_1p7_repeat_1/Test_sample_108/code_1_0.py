import math

def solve_moment_capacity():
    """
    Calculates the moment at collapse for the given non-standard section.
    """
    # 1. Inputs and Constants based on the problem description
    fck = 30  # MPa
    fyk = 500  # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    E_s = 200000  # MPa (Modulus of Elasticity for steel)
    eps_cu3 = 0.0035  # Ultimate compressive strain in concrete

    # Geometry (all units in mm)
    h_trap = 300  # Height of the top trapezoidal part
    b_top = 100   # Width at the very top
    b_bottom = 400 # Width at the base of the trapezoid and the section bottom
    
    # Reinforcement details
    d_top = 50                 # Depth of top reinforcement
    d_mid = 50 + 210           # Depth of middle reinforcement
    d_bot = 50 + 210 + 90      # Depth of bottom reinforcement
    d_bar = 20                 # Diameter of reinforcement bars

    # 2. Derived Material Properties and Areas
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    eps_yd = fyd / E_s
    A_bar = math.pi * (d_bar / 2)**2
    As_top = 2 * A_bar  # 2H20
    As_mid = 2 * A_bar  # 2H20
    As_bot = 3 * A_bar  # 3H20

    print("--- Design Parameters ---")
    print(f"fcd (Design strength of concrete) = {fcd:.2f} MPa")
    print(f"fyd (Design strength of steel) = {fyd:.2f} MPa")
    print(f"Area of one H20 bar = {A_bar:.2f} mm^2")
    print("-" * 25)

    # 3. Iteratively find the Plastic Neutral Axis (x)
    def get_force_imbalance(x):
        # Calculates F_compression - F_tension for a given x
        if x <= 0: return -1

        s = 0.8 * x  # Depth of rectangular stress block

        # Concrete force calculation (for x within the top trapezoidal part)
        # Width b(y) = 100 + y for 0 <= y <= 300
        if s <= h_trap:
            Ac = b_top * s + s**2 / 2
            Fc = fcd * Ac
        else: # PNA is in the lower rectangular part (unlikely for this problem but included for completeness)
            Ac_trap = b_top * h_trap + h_trap**2 / 2
            Fc_trap = fcd * Ac_trap
            Fc_rect = fcd * b_bottom * (s - h_trap)
            Fc = Fc_trap + Fc_rect

        # Steel forces
        # Top steel (compression)
        eps_s_top = eps_cu3 * (x - d_top) / x
        sigma_s_top = -fyd if abs(eps_s_top) > eps_yd else E_s * eps_s_top
        F_s_top = As_top * sigma_s_top  # Will be negative

        # Middle steel (tension)
        eps_s_mid = eps_cu3 * (d_mid - x) / x
        sigma_s_mid = fyd if eps_s_mid > eps_yd else (E_s * eps_s_mid if eps_s_mid > 0 else 0)
        F_s_mid = As_mid * sigma_s_mid

        # Bottom steel (tension)
        eps_s_bot = eps_cu3 * (d_bot - x) / x
        sigma_s_bot = fyd if eps_s_bot > eps_yd else (E_s * eps_s_bot if eps_s_bot > 0 else 0)
        F_s_bot = As_bot * sigma_s_bot
        
        F_comp_total = Fc - F_s_top # F_s_top is negative
        F_tens_total = F_s_mid + F_s_bot
        
        return F_comp_total - F_tens_total

    # Bisection search for x
    x_low, x_high = 1.0, h_trap / 0.8
    for _ in range(100):
        x_guess = (x_low + x_high) / 2
        imbalance = get_force_imbalance(x_guess)
        if imbalance > 0: # Compression too high, x is too deep
            x_high = x_guess
        else: # Tension too high, x is too shallow
            x_low = x_guess
    x = (x_low + x_high) / 2

    print(f"--- Neutral Axis Calculation ---")
    print(f"Found plastic neutral axis depth x = {x:.2f} mm")
    print("-" * 25)

    # 4. Calculate Collapse Moment (M_rd)
    s = 0.8 * x

    # Final forces and their locations
    # Concrete force and centroid
    Ac = b_top * s + s**2 / 2
    Fc = fcd * Ac
    moment_of_area = b_top/2 * s**2 + s**3 / 3
    yc = moment_of_area / Ac

    # Steel forces
    eps_s_top = eps_cu3 * (x - d_top) / x
    sigma_s_top = -fyd if abs(eps_s_top) > eps_yd else E_s * eps_s_top
    F_s_top = As_top * sigma_s_top
    
    eps_s_mid = eps_cu3 * (d_mid - x) / x
    sigma_s_mid = fyd if eps_s_mid > eps_yd else E_s * eps_s_mid
    F_s_mid = As_mid * sigma_s_mid

    eps_s_bot = eps_cu3 * (d_bot - x) / x
    sigma_s_bot = fyd if eps_s_bot > eps_yd else E_s * eps_s_bot
    F_s_bot = As_bot * sigma_s_bot
    
    F_comp_res = Fc - F_s_top
    F_tens_res = F_s_mid + F_s_bot
    
    # Location of resultant forces from top fiber
    d_comp = (Fc * yc - F_s_top * d_top) / F_comp_res
    d_tens = (F_s_mid * d_mid + F_s_bot * d_bot) / F_tens_res
    
    # Internal lever arm and Moment
    z = d_tens - d_comp
    M_rd = F_tens_res * z
    M_rd_kNm = M_rd / 1e6
    
    print("--- Collapse Moment Calculation ---")
    print("\nIndividual Forces (1 kN = 1000 N):")
    print(f"  Fc (Concrete compression)       = {Fc/1000:.2f} kN at yc = {yc:.2f} mm")
    print(f"  Fs_top (Compression steel)      = {abs(F_s_top)/1000:.2f} kN at d = {d_top:.1f} mm")
    print(f"  Fs_mid (Tension steel)          = {F_s_mid/1000:.2f} kN at d = {d_mid:.1f} mm")
    print(f"  Fs_bot (Tension steel)          = {F_s_bot/1000:.2f} kN at d = {d_bot:.1f} mm")
    
    print("\nResultant Forces and Locations:")
    print(f"  Resultant Tension Force F_tens  = {F_tens_res/1000:.2f} kN")
    print(f"  Resultant Compression Force F_comp = {F_comp_res/1000:.2f} kN")
    print(f"  Location of F_tens, d_tens      = ({F_s_mid/1000:.2f} * {d_mid:.1f} + {F_s_bot/1000:.2f} * {d_bot:.1f}) / {F_tens_res/1000:.2f} = {d_tens:.2f} mm")
    print(f"  Location of F_comp, d_comp      = ({Fc/1000:.2f} * {yc:.2f} + {abs(F_s_top)/1000:.2f} * {d_top:.1f}) / {F_comp_res/1000:.2f} = {d_comp:.2f} mm")

    print("\nFinal Moment Calculation:")
    print(f"  Internal Lever Arm z = d_tens - d_comp = {d_tens:.2f} - {d_comp:.2f} = {z:.2f} mm")
    print(f"  Moment M_rd = F_tens * z = {F_tens_res/1000:.2f} kN * {z:.2f} mm = {M_rd/1e6:.2f} kNm")

    return M_rd_kNm

final_moment = solve_moment_capacity()
print(f"\nFinal Answer: The moment at collapse is {final_moment:.2f} kNm.")
print(f"<<<{final_moment:.1f}>>>")
