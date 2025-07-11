import math

def solve_moment_capacity():
    """
    Calculates the moment at collapse for the given non-standard concrete section.
    """
    # 1. Define Material Properties and Safety Factors (Eurocode 2)
    fck = 30.0  # MPa
    fyk = 500.0  # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    Es = 200000.0  # MPa (Young's Modulus for steel)
    lambda_val = 0.8 # Factor for depth of rectangular stress block
    epsilon_cu = 0.0035 # Ultimate concrete strain

    # Design strengths
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    epsilon_yd = fyd / Es

    print("--- Material Design Properties ---")
    print(f"Design compressive strength of concrete (fcd): {fcd:.2f} MPa")
    print(f"Design yield strength of steel (fyd): {fyd:.2f} MPa")
    print(f"Yield strain of steel (ε_yd): {epsilon_yd:.5f}\n")

    # 2. Define Section Geometry and Reinforcement
    bar_diameter = 20.0  # mm
    as_bar = math.pi * (bar_diameter / 2)**2

    # Reinforcement layers (area and depth from top edge)
    # Compression steel
    As_comp = 2 * as_bar
    d_prime = 50.0  # mm
    # Tension steel layer 1
    As_tens1 = 2 * as_bar
    d1 = 260.0  # mm
    # Tension steel layer 2
    As_tens2 = 3 * as_bar
    d2 = 350.0  # mm

    print("--- Reinforcement Details ---")
    print(f"Area of compression steel (As'): {As_comp:.2f} mm^2 at d' = {d_prime} mm")
    print(f"Area of tension steel layer 1 (As1): {As_tens1:.2f} mm^2 at d1 = {d1} mm")
    print(f"Area of tension steel layer 2 (As2): {As_tens2:.2f} mm^2 at d2 = {d2} mm\n")

    # 3. & 4. Determine the Neutral Axis Depth (x)
    # The equilibrium equation C=T leads to a cubic equation for x:
    # A*x^3 + B*x^2 + C*x + D = 0
    # where forces are calculated based on the assumption that compression steel and
    # the bottom layer of tension steel yield, but the middle layer does not.
    
    # Coefficients of the cubic equation: 5.44x^3 + 1360x^2 + 292547x - 114353967 = 0
    def equation(x):
        Fsc_net = As_comp * (fyd - fcd)
        Fst2 = As_tens2 * fyd
        # Stress in non-yielded steel
        sigma_st1 = Es * epsilon_cu * (d1 - x) / x
        Fst1 = As_tens1 * sigma_st1
        
        s = lambda_val * x
        # Area of concrete in compression (trapezoid)
        Ac = 100 * s + 0.5 * s**2
        Cc = fcd * Ac
        
        return Cc + Fsc_net - Fst1 - Fst2

    # Solve for x using bisection method
    low, high = 100, 200
    x = (low + high) / 2
    for _ in range(30): # 30 iterations for high precision
        if equation(x) * equation(high) < 0:
            low = x
        else:
            high = x
        x = (low + high) / 2
        
    print("--- Neutral Axis Calculation ---")
    print(f"Solving for neutral axis depth 'x' by balancing forces C=T.")
    print(f"Found neutral axis depth (x): {x:.2f} mm\n")

    # Verify strain assumptions
    epsilon_sc = epsilon_cu * (x - d_prime) / x
    epsilon_st1 = epsilon_cu * (d1 - x) / x
    epsilon_st2 = epsilon_cu * (d2 - x) / x
    
    print("--- Strain Verification ---")
    print(f"Strain in compression steel (ε_sc): {epsilon_sc:.5f} (> ε_yd={epsilon_yd:.5f}, so it yields)")
    print(f"Strain in tension steel layer 1 (ε_st1): {epsilon_st1:.5f} (< ε_yd={epsilon_yd:.5f}, so it does NOT yield)")
    print(f"Strain in tension steel layer 2 (ε_st2): {epsilon_st2:.5f} (> ε_yd={epsilon_yd:.5f}, so it yields)\n")

    # 5. Calculate Internal Forces and Lever Arms
    s = lambda_val * x
    
    # Concrete force and centroid
    Ac = 100 * s + 0.5 * s**2
    Cc = fcd * Ac
    moment_Ac = 50 * s**2 + (s**3) / 3
    yc = moment_Ac / Ac

    # Steel forces
    Fsc_net = As_comp * (fyd - fcd) # Net force, accounting for displaced concrete
    sigma_st1 = Es * epsilon_st1
    Fst1 = As_tens1 * sigma_st1
    Fst2 = As_tens2 * fyd

    print("--- Internal Forces and Centroid ---")
    print(f"Concrete compressive force (Cc): {Cc/1000:.2f} kN")
    print(f"Centroid of concrete force (yc): {yc:.2f} mm from top")
    print(f"Net compression steel force (Fsc_net): {Fsc_net/1000:.2f} kN")
    print(f"Tension steel force layer 1 (Fst1): {Fst1/1000:.2f} kN")
    print(f"Tension steel force layer 2 (Fst2): {Fst2/1000:.2f} kN\n")

    # 6. Calculate the Moment at Collapse (Mu)
    # Sum of moments about the top edge of the section
    moment_tension = Fst1 * d1 + Fst2 * d2
    moment_compression = Fsc_net * d_prime + Cc * yc
    
    Mu_Nmm = moment_tension - moment_compression
    Mu_kNm = Mu_Nmm / 1e6

    print("--- Moment Capacity Calculation ---")
    print("Mu = (Fst1 * d1 + Fst2 * d2) - (Fsc_net * d' + Cc * yc)")
    print(f"Mu = ({Fst1/1000:.2f} kN * {d1} mm + {Fst2/1000:.2f} kN * {d2} mm) - ({Fsc_net/1000:.2f} kN * {d_prime} mm + {Cc/1000:.2f} kN * {yc:.2f} mm)")
    print(f"Mu = ({moment_tension/1e6:.2f}) - ({moment_compression/1e6:.2f}) kNm")
    print(f"Mu = {Mu_kNm:.2f} kNm\n")
    
    return Mu_kNm

# Run the calculation and print the final answer
final_moment = solve_moment_capacity()
print(f"The moment at collapse for the section is {final_moment:.1f} kNm.")
print(f"<<<{final_moment:.1f}>>>")
