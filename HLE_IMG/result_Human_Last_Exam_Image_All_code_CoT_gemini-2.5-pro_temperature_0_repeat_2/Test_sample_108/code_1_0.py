import math

def calculate_moment_at_collapse():
    """
    Calculates the moment at collapse for the given non-standard concrete section.
    """
    # --- 1. Initial parameters and constants ---
    fck = 30.0  # MPa, characteristic compressive strength of concrete
    fyk = 500.0  # MPa, characteristic yield strength of steel
    gamma_c = 1.5  # Partial safety factor for concrete
    gamma_s = 1.15  # Partial safety factor for steel
    alpha_cc = 0.85  # Coefficient for long-term effects on concrete strength
    E_s = 200000.0  # MPa, Young's modulus for steel
    eps_cu = 0.0035  # Ultimate compressive strain in concrete
    lambda_val = 0.8  # Factor for depth of rectangular stress block
    eta = 1.0 # Factor for effective strength in stress block

    # --- 2. Derived material properties ---
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    eps_yd = fyd / E_s
    d_bar = 20.0 # mm, bar diameter
    A_s20 = math.pi * (d_bar / 2)**2

    print("--- Material Properties ---")
    print(f"Design compressive strength of concrete, fcd = {alpha_cc} * {fck} / {gamma_c} = {fcd:.2f} MPa")
    print(f"Design yield strength of steel, fyd = {fyk} / {gamma_s} = {fyd:.2f} MPa")
    print(f"Area of one H20 bar, A_s20 = pi * (20/2)^2 = {A_s20:.2f} mm^2\n")

    # --- 3. Reinforcement details ---
    n_top = 2
    A_s_top = n_top * A_s20
    d_top = 50.0  # mm

    n_mid = 2
    A_s_mid = n_mid * A_s20
    d_mid = 50.0 + 210.0

    n_bot = 3
    A_s_bot = n_bot * A_s20
    d_bot = 50.0 + 210.0 + 90.0

    A_st_total = A_s_mid + A_s_bot

    print("--- Reinforcement Details ---")
    print(f"Top compression steel: {n_top}H20, Area A_s_top = {A_s_top:.2f} mm^2 at d_top = {d_top:.1f} mm")
    print(f"Middle tension steel: {n_mid}H20, Area A_s_mid = {A_s_mid:.2f} mm^2 at d_mid = {d_mid:.1f} mm")
    print(f"Bottom tension steel: {n_bot}H20, Area A_s_bot = {A_s_bot:.2f} mm^2 at d_bot = {d_bot:.1f} mm\n")

    # --- 4. Force Equilibrium to find Neutral Axis (x) ---
    print("--- Finding Neutral Axis (x) via Force Equilibrium ---")
    print("Equilibrium Equation: F_concrete + F_compression_steel = F_tension_steel")
    
    # Total tension force (assuming all tension steel yields)
    F_st = A_st_total * fyd
    print(f"Total Tension Force (assuming yield), F_st = {A_st_total:.2f} * {fyd:.2f} = {F_st:.0f} N")

    # Iterative calculation to find x
    x = d_bot / 2  # Initial guess
    for i in range(10):
        # Strain and stress in compression steel
        eps_sc = eps_cu * (x - d_top) / x
        sigma_sc = min(E_s * eps_sc, fyd)
        
        # Force in compression steel (accounting for displaced concrete)
        F_sc = (sigma_sc - eta * fcd) * A_s_top
        if F_sc < 0: F_sc = 0

        # Required concrete compression force
        F_c_req = F_st - F_sc
        
        # Required area of concrete compression block
        A_c_req = F_c_req / (eta * fcd)
        
        # From A_c, find depth of stress block 's' for the trapezoid
        # A_c = 100*s + 0.5*s^2  => 0.5*s^2 + 100*s - A_c = 0
        a_q, b_q, c_q = 0.5, 100.0, -A_c_req
        s = (-b_q + math.sqrt(b_q**2 - 4 * a_q * c_q)) / (2 * a_q)
        
        x_new = s / lambda_val
        if abs(x_new - x) < 0.01:
            x = x_new
            break
        x = x_new

    print(f"Converged Neutral Axis Depth, x = {x:.2f} mm\n")

    # --- 5. Verification of Assumptions ---
    eps_st_bot = eps_cu * (d_bot - x) / x
    print("--- Verification of Assumptions ---")
    print(f"Strain in bottom tension steel = {eps_cu:.4f} * ({d_bot:.1f} - {x:.2f}) / {x:.2f} = {eps_st_bot:.5f}")
    print(f"Yield strain of steel = {eps_yd:.5f}")
    if eps_st_bot >= eps_yd:
        print("-> Tension steel has yielded. Assumption is valid.\n")
    else:
        print("-> Tension steel has NOT yielded. Analysis is invalid.\n")

    # --- 6. Moment Calculation ---
    print("--- Moment at Collapse Calculation ---")
    print("Taking moments of internal forces about the top fiber.")
    
    # Final forces based on converged x
    s = lambda_val * x
    eps_sc = eps_cu * (x - d_top) / x
    sigma_sc = min(E_s * eps_sc, fyd)
    F_sc = (sigma_sc - eta * fcd) * A_s_top
    if F_sc < 0: F_sc = 0
    
    A_c = 100*s + 0.5*s**2
    F_c = eta * fcd * A_c
    F_st_mid = A_s_mid * fyd
    F_st_bot = A_s_bot * fyd

    # Centroid of concrete compression area (z_c) from top fiber
    moment_of_area = 50 * s**2 + (s**3) / 3
    z_c = moment_of_area / A_c

    # Moment calculation
    M_tension = F_st_mid * d_mid + F_st_bot * d_bot
    M_compression = F_c * z_c + F_sc * d_top
    M_Rd_Nmm = M_tension - M_compression
    M_Rd_kNm = M_Rd_Nmm / 1e6

    print("\nMoment Equation: M_Rd = (F_st_mid * d_mid + F_st_bot * d_bot) - (F_c * z_c + F_sc * d_top)")
    print(f"M_Rd = ({F_st_mid:.0f} N * {d_mid:.1f} mm + {F_st_bot:.0f} N * {d_bot:.1f} mm) - ({F_c:.0f} N * {z_c:.2f} mm + {F_sc:.0f} N * {d_top:.1f} mm)")
    print(f"M_Rd = ({M_tension:.0f}) - ({M_compression:.0f})")
    print(f"M_Rd = {M_Rd_Nmm:.0f} Nmm")
    
    print("\n--- Final Answer ---")
    print(f"The moment at collapse for the section is {M_Rd_kNm:.2f} kNm.")
    return M_Rd_kNm

if __name__ == '__main__':
    moment = calculate_moment_at_collapse()
    # The final answer is extracted from the last print statement.
    # For the final submission format:
    # print(f"<<<{moment:.2f}>>>")