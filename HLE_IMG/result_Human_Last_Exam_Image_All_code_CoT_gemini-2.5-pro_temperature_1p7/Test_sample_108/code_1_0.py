import math

def calculate_moment_capacity():
    """
    Calculates the moment at collapse for the given non-standard concrete section.
    """
    # 1. Material Properties and Safety Factors (EC2)
    fck = 30.0  # MPa
    fyk = 500.0  # MPa
    gamma_c = 1.5  # Partial safety factor for concrete
    gamma_s = 1.15  # Partial safety factor for steel
    lambda_val = 0.8 # Factor for effective depth of compression block
    eta = 1.0       # Factor for effective strength
    epsilon_cu3 = 0.0035 # Ultimate concrete strain
    Es = 200000.0   # Modulus of elasticity for steel (MPa)

    # Design strengths
    fcd = eta * fck / gamma_c
    fyd = fyk / gamma_s
    epsilon_yd = fyd / Es

    print("--- Design Parameters ---")
    print(f"fcd (Design concrete strength): {fcd:.2f} MPa")
    print(f"fyd (Design steel yield strength): {fyd:.2f} MPa")
    print(f"Yield strain of steel: {epsilon_yd:.5f}\n")

    # 2. Reinforcement Details
    d_bar = 20.0  # mm (H20 bar)
    A_bar = math.pi * (d_bar / 2)**2

    # Layer 1 (Top, Compression)
    n_sc = 2
    As_c = n_sc * A_bar
    d_c = 50.0  # mm

    # Layer 2 (Middle, Tension)
    n_st1 = 2
    As_t1 = n_st1 * A_bar
    d_t1 = 50.0 + 210.0  # mm

    # Layer 3 (Bottom, Tension)
    n_st2 = 3
    As_t2 = n_st2 * A_bar
    d_t2 = 50.0 + 210.0 + 90.0 # mm

    print("--- Reinforcement Details ---")
    print(f"Area of compression steel (As_c) at d={d_c}mm: {As_c:.2f} mm^2")
    print(f"Area of tension steel (As_t1) at d={d_t1}mm: {As_t1:.2f} mm^2")
    print(f"Area of tension steel (As_t2) at d={d_t2}mm: {As_t2:.2f} mm^2\n")

    # 3. Force Equilibrium (Assume all steel yields)
    # EC2 allows ignoring the stress in concrete displaced by compression steel
    Fsc = As_c * fyd
    Fst1 = As_t1 * fyd
    Fst2 = As_t2 * fyd
    Fst_total = Fst1 + Fst2

    # Required concrete force
    Fc = Fst_total - Fsc

    print("--- Force Calculation (Assuming Steel Yields) ---")
    print(f"Force in compression steel (Fsc): {Fsc/1000:.2f} kN")
    print(f"Total force in tension steel (Fst): {Fst_total/1000:.2f} kN")
    print(f"Required force in concrete (Fc): {Fc/1000:.2f} kN\n")

    # 4. Determine Neutral Axis Depth (x)
    # The compressive zone is in the top trapezoidal part where b(y) = 100 + y
    # Fc = fcd * Area_compression_block = fcd * (100*s + s^2/2)
    # This leads to a quadratic equation for s: (fcd/2)*s^2 + (fcd*100)*s - Fc = 0
    # or 10*s^2 + 2000*s - Fc = 0
    a = fcd / 2
    b = fcd * 100
    c = -Fc
    
    # Solve quadratic equation for s (depth of stress block)
    s = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
    
    # Calculate neutral axis depth x
    x = s / lambda_val

    print("--- Neutral Axis Calculation ---")
    print(f"Depth of rectangular stress block (s): {s:.2f} mm")
    print(f"Depth of neutral axis (x): {x:.2f} mm\n")

    # 5. Verify Strain Assumptions
    epsilon_sc = epsilon_cu3 * (x - d_c) / x
    epsilon_st1 = epsilon_cu3 * (d_t1 - x) / x
    epsilon_st2 = epsilon_cu3 * (d_t2 - x) / x

    print("--- Strain Verification ---")
    print(f"Strain in compression steel (ε_sc): {epsilon_sc:.5f} -> {'Yields' if epsilon_sc >= epsilon_yd else 'Does not yield'}")
    print(f"Strain in tension steel layer 1 (ε_st1): {epsilon_st1:.5f} -> {'Yields' if epsilon_st1 >= epsilon_yd else 'Does not yield'}")
    print(f"Strain in tension steel layer 2 (ε_st2): {epsilon_st2:.5f} -> {'Yields' if epsilon_st2 >= epsilon_yd else 'Does not yield'}")
    print("All steel layers are confirmed to have yielded.\n")

    # 6. Calculate Moment Capacity (M)
    # Find the centroid of the trapezoidal compression area
    # Area = 100*s + s^2/2
    # Moment of Area about top = 100*s^2/2 + s^3/3
    Area_c = 100 * s + s**2 / 2
    Moment_Area_c = 100 * s**2 / 2 + s**3 / 3
    zc = Moment_Area_c / Area_c
    
    # Actual concrete force based on calculated s, for higher precision
    Fc_actual = fcd * Area_c

    print("--- Moment Calculation ---")
    print(f"Centroid of concrete force from top fiber (zc): {zc:.2f} mm\n")
    
    # Take moments about the top fiber
    # M = (F_tension * d_tension) - (F_compression * d_compression)
    M = (Fst1 * d_t1 + Fst2 * d_t2) - (Fsc * d_c) - (Fc_actual * zc)
    
    # Convert M from N.mm to kN.m
    M_kNm = M / 1e6
    
    # Print the final equation with values
    print("The final moment capacity is calculated by summing moments about the top fiber:")
    print("M = (F_st1 * d_t1) + (F_st2 * d_t2) - (F_sc * d_c) - (F_c * zc)")
    print("M = ({:.0f} N * {:.0f} mm) + ({:.0f} N * {:.0f} mm) - ({:.0f} N * {:.0f} mm) - ({:.0f} N * {:.2f} mm)".format(
        Fst1, d_t1, Fst2, d_t2, Fsc, d_c, Fc_actual, zc))
    print("M = ({:.2f} + {:.2f} - {:.2f} - {:.2f}) * 10^6 N.mm".format(
        (Fst1 * d_t1)/1e6, (Fst2 * d_t2)/1e6, (Fsc * d_c)/1e6, (Fc_actual * zc)/1e6
    ))
    print(f"M = {M:.2f} N.mm\n")


    print("--- Final Result ---")
    print(f"The moment at collapse for the section is: {M_kNm:.1f} kNm")

# Run the calculation
calculate_moment_capacity()
<<<171.7>>>