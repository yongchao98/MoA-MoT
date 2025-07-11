import math

def solve_moment_capacity():
    """
    Calculates the moment at collapse for the given non-standard concrete section.
    """
    # 1. Define Material Properties and Constants (Eurocode 2)
    fck = 30.0  # MPa, characteristic compressive cylinder strength of concrete
    fyk = 500.0  # MPa, characteristic yield strength of reinforcement
    gamma_c = 1.5  # Partial safety factor for concrete
    gamma_s = 1.15 # Partial safety factor for steel
    alpha_cc = 0.85 # Coefficient for long-term effects
    lambda_ = 0.8  # Factor for effective height of compression block
    eta = 1.0     # Factor for effective strength
    Es = 200000.0 # MPa, Modulus of elasticity of steel
    ecu = 0.0035  # Ultimate concrete strain in compression

    # Calculate design strengths
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    eyd = fyd / Es # Yield strain of steel

    print("--- Material Design Properties ---")
    print(f"fcd (Design compressive strength of concrete) = {fcd:.2f} MPa")
    print(f"fyd (Design yield strength of steel) = {fyd:.2f} MPa")
    print(f"eyd (Design yield strain of steel) = {eyd:.5f}\n")

    # 2. Define Geometry and Reinforcement (units in mm, mm^2)
    # Steel Layers from top edge
    # Layer 1: 2H20 (Compression)
    d1 = 50.0
    As1 = 2 * math.pi * (20/2)**2
    # Layer 2: 2H20 (Tension)
    d2 = 50.0 + 210.0
    As2 = 2 * math.pi * (20/2)**2
    # Layer 3: 3H20 (Tension)
    d3 = 50.0 + 210.0 + 90.0
    As3 = 3 * math.pi * (20/2)**2

    print("--- Reinforcement Details ---")
    print(f"As1 (top, compression) at d1={d1} mm: {As1:.2f} mm^2")
    print(f"As2 (middle, tension) at d2={d2} mm: {As2:.2f} mm^2")
    print(f"As3 (bottom, tension) at d3={d3} mm: {As3:.2f} mm^2\n")
    
    # 3. Iterative Solution to find Neutral Axis depth 'x'
    
    def get_forces_and_residual(x):
        # Calculates internal forces and their residual for a given neutral axis depth 'x'.

        # Assume x > d1
        if x <= d1:
            return 1e9 # Physically unreasonable, return large residual

        # --- Calculate Stresses and Forces in Steel ---
        # Compression steel (Layer 1)
        esc = ecu * (x - d1) / x
        # Stress in compression steel (capped at fyd)
        sigma_sc = min(Es * esc, fyd)
        # Force in compression steel
        Fsc = As1 * sigma_sc

        # Tension steel (Layer 2)
        if x < d2:
            est2 = ecu * (d2 - x) / x
            sigma_st2 = min(Es * est2, fyd)
        else: # Layer is in compression
            est2 = ecu * (x - d2) / x
            sigma_st2 = -min(Es * est2, fyd) # Negative for tension force calculation
        Fst2 = As2 * sigma_st2

        # Tension steel (Layer 3)
        est3 = ecu * (d3 - x) / x
        # Stress in tension steel (capped at fyd)
        sigma_st3 = min(Es * est3, fyd)
        Fst3 = As3 * sigma_st3
        
        # --- Calculate Force in Concrete Compression Block ---
        s = lambda_ * x  # Effective depth of compression block
        
        # Geometry of the section: Top trapezoid (h=300, w=100->400), Bottom rect (h=100, w=400)
        # Width at a depth y from the top (for y <= 300) is b(y) = 100 + y
        if s <= 300: # Compression block is a trapezoid
            width_top = 100.0
            width_bottom = 100.0 + s
            Ac = (width_top + width_bottom) / 2.0 * s
        else: # Compression block covers top trapezoid and part of bottom rectangle
            Area_trap_full = (100.0 + 400.0) / 2.0 * 300.0
            Area_rect = 400.0 * (s - 300.0)
            Ac = Area_trap_full + Area_rect
            
        Fc = eta * fcd * Ac

        # --- Total Forces and Residual ---
        Total_Compression = Fc + Fsc
        Total_Tension = Fst2 + Fst3
        residual = Total_Compression - Total_Tension
        
        forces = {
            'Fc': Fc, 'Fsc': Fsc, 'Fst2': Fst2, 'Fst3': Fst3,
            'sigma_sc': sigma_sc, 'sigma_st2': sigma_st2, 'sigma_st3': sigma_st3
        }
        return residual, forces

    # Bisection method to find x
    x_low = d1 # Lower bound for x
    x_high = d3 # Upper bound for x
    
    # Check bounds
    res_low, _ = get_forces_and_residual(x_low)
    res_high, _ = get_forces_and_residual(x_high)
    if res_low * res_high > 0:
        print("Error: The root is not bracketed. Check bounds or function.")
        return

    for i in range(100): # 100 iterations for high precision
        x_mid = (x_low + x_high) / 2.0
        residual, _ = get_forces_and_residual(x_mid)
        if residual > 0: # C > T, so x is too large
            x_high = x_mid
        else: # C < T, so x is too small
            x_low = x_mid
        if abs(residual) < 1.0: # Convergence tolerance (1 N)
            break
            
    x = (x_low + x_high) / 2.0
    residual, forces = get_forces_and_residual(x)

    print(f"--- Equilibrium Analysis ---")
    print(f"Final Neutral Axis Depth, x = {x:.2f} mm")
    print(f"C (Compression) = {forces['Fc']/1000:.2f} + {forces['Fsc']/1000:.2f} = {(forces['Fc']+forces['Fsc'])/1000:.2f} kN")
    print(f"T (Tension)     = {forces['Fst2']/1000:.2f} + {forces['Fst3']/1000:.2f} = {(forces['Fst2']+forces['Fst3'])/1000:.2f} kN")
    print(f"Residual C-T = {residual:.2f} N\n")

    # 4. Calculate Moment at Collapse (M_rd)
    # Take moments about the top compressive fiber (y=0)
    
    # Moment from steel forces
    M_sc = forces['Fsc'] * d1
    M_st2 = forces['Fst2'] * d2
    M_st3 = forces['Fst3'] * d3

    # Moment from concrete force
    s = lambda_ * x
    Fc = forces['Fc']
    if s <= 300: # Centroid of trapezoidal compression block
        a = 100.0  # Top width
        b = 100.0 + s  # Bottom width
        h = s
        # Centroid distance from top edge (side 'a')
        zc = (h / 3.0) * (a + 2*b) / (a + b)
    else: # Composite centroid (trapezoid + rectangle)
        Area_trap_full = (100.0 + 400.0) / 2.0 * 300.0
        a_trap = 100; b_trap = 400; h_trap = 300
        zc_trap = (h_trap / 3.0) * (a_trap + 2*b_trap) / (a_trap + b_trap)
        
        Area_rect = 400.0 * (s - 300.0)
        zc_rect = 300.0 + (s - 300.0)/2.0
        
        zc = (Area_trap_full * zc_trap + Area_rect * zc_rect) / (Area_trap_full + Area_rect)

    Mc = Fc * zc

    # Total moment
    M_rd_Nmm = M_st2 + M_st3 - M_sc - Mc
    M_rd_kNm = M_rd_Nmm / 1e6

    print("--- Moment Calculation (about top fiber) ---")
    print(f"Lever arm of Fc, zc = {zc:.2f} mm")
    print("Final Moment Equation: M_rd = F_st2*d2 + F_st3*d3 - F_sc*d1 - F_c*zc")
    print(f"M_rd = {forces['Fst2']/1000:.2f}kN * {d2/1000:.3f}m + {forces['Fst3']/1000:.2f}kN * {d3/1000:.3f}m - {forces['Fsc']/1000:.2f}kN * {d1/1000:.3f}m - {forces['Fc']/1000:.2f}kN * {zc/1000:.3f}m")
    print(f"M_rd = {(M_st2)/1e6:.2f} + {(M_st3)/1e6:.2f} - {(M_sc)/1e6:.2f} - {(Mc)/1e6:.2f} kNm")
    print(f"Moment at collapse, M_rd = {M_rd_kNm:.2f} kNm")
    
    return M_rd_kNm

# Execute the function and get the final answer
moment_capacity = solve_moment_capacity()
print(f"\nFinal Answer: The moment at collapse is {moment_capacity:.2f} kNm.")
print(f"<<<{moment_capacity:.1f}>>>")
