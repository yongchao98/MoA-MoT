import math

def solve_footing_resistance():
    """
    Calculates the required ULS design bearing resistance for a pad footing
    by iteratively finding the minimum required footing size.
    """
    
    # --- 1. Initial Parameters from the problem statement ---
    Gk = 1000.0  # kN, Characteristic Permanent Load
    Qk_v = 1500.0 # kN, Characteristic Variable Vertical Load
    Qk_h = 300.0  # kN, Characteristic Variable Horizontal Load
    h_load_app = 2.0 # m, height of horizontal load application above ground
    Df = 0.75 # m, footing depth
    gamma_soil = 20.0 # kN/m^3, soil unit weight
    phi_k_deg = 35.0 # degrees, characteristic soil friction angle
    c_prime_k = 0.0 # kPa, characteristic effective cohesion

    # --- 2. ULS Combination 1 Design Values (EC7, DA1-1) ---
    # Partial factors for actions (Set A1) and soil properties (Set M1)
    gamma_G = 1.35
    gamma_Q = 1.50
    gamma_phi = 1.0

    # Design loads and moments
    Vd = gamma_G * Gk + gamma_Q * Qk_v
    Hd = gamma_Q * Qk_h
    h_total = h_load_app + Df
    Md = Hd * h_total
    e = Md / Vd # Load Eccentricity

    # --- 3. Design Soil Parameters ---
    phi_k_rad = math.radians(phi_k_deg)
    phi_d_rad = math.atan(math.tan(phi_k_rad) / gamma_phi)
    c_prime_d = c_prime_k / 1.0 # using gamma_c=1.0
    q_prime = gamma_soil * Df # effective overburden pressure

    # --- 4. Bearing Capacity Factors (from EC7 Annex D) ---
    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
    Nc = (Nq - 1) / math.tan(phi_d_rad)
    Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)

    # --- 5. Iterative Search for Minimum Footing Width B ---
    B_min_uplift = 6 * e # Theoretical minimum B to avoid uplift
    B = B_min_uplift
    step = 0.001 # Search step for B in meters

    while True:
        B += step
        L = B # Square footing B x L

        # Effective dimensions considering eccentricity
        B_prime = B - 2 * e
        if B_prime <= 0:
            continue # This B is too small, load is outside the footing
        
        L_prime = L
        A_prime = B_prime * L_prime
        
        # ULS Applied Pressure on effective area
        q_applied = Vd / A_prime

        # Shape factors (s)
        ratio_BL = B_prime / L_prime
        s_q = 1 + ratio_BL * math.sin(phi_d_rad)
        s_gamma = 1 - 0.3 * ratio_BL
        
        # Load inclination factors (i) - assuming H is parallel to B'
        m = (2 + ratio_BL + 1) / (2 * ratio_BL + 1)
        i_q = (1 - Hd / Vd)**m
        i_gamma = (1 - Hd / Vd)**(m + 1)
        
        # ULS Design Bearing Resistance (per unit effective area)
        # Note: The term with c_prime_d is zero since c' = 0
        term_q = q_prime * Nq * s_q * i_q
        term_gamma = 0.5 * gamma_soil * B_prime * Ngamma * s_gamma * i_gamma
        q_resistance = term_q + term_gamma
        
        # Check if resistance is sufficient
        if q_resistance >= q_applied:
            # Found the minimum required B
            required_resistance = q_applied
            break

    # --- 6. Output Results ---
    print(f"The required ULS design bearing resistance is the pressure exerted by the footing when sized at the limit of stability.\n")
    print(f"--- Design Calculation Summary ---")
    print(f"1. ULS Design Loads & Eccentricity:")
    print(f"   Vd = {gamma_G} * {Gk} + {gamma_Q} * {Qk_v} = {Vd:.1f} kN")
    print(f"   Hd = {gamma_Q} * {Qk_h} = {Hd:.1f} kN")
    print(f"   e = ({Hd:.1f} kN * ({h_load_app} + {Df}) m) / {Vd:.1f} kN = {e:.4f} m")

    print(f"\n2. Sizing Results (at minimum required size):")
    print(f"   An iterative search finds the minimum footing width B ≈ {B:.3f} m.")
    print(f"   For this size:")
    print(f"   Effective width B' = {B:.3f} m - 2 * {e:.4f} m = {B_prime:.4f} m")
    print(f"   Effective area A' = {B_prime:.4f} m * {B:.3f} m = {A_prime:.4f} m^2")
    
    print(f"\n3. Applied Pressure vs. Soil Resistance (at B={B:.3f} m):")
    print(f"   Applied Pressure (Demand):")
    print(f"   q_applied = Vd / A' = {Vd:.1f} kN / {A_prime:.4f} m^2 = {q_applied:.1f} kN/m^2")

    print(f"\n   Soil Resistance (Capacity):")
    print(f"   q_resistance = (q' * Nq * s_q * i_q) + (0.5 * γ' * B' * Nγ * s_γ * i_γ)")
    print(f"   q'_term      = {q_prime:.1f} * {Nq:.2f} * {s_q:.3f} * {i_q:.3f} = {term_q:.1f} kN/m^2")
    print(f"   γ'_term      = 0.5 * {gamma_soil:.1f} * {B_prime:.4f} * {Ngamma:.2f} * {s_gamma:.3f} * {i_gamma:.3f} = {term_gamma:.1f} kN/m^2")
    print(f"   q_resistance = {term_q:.1f} + {term_gamma:.1f} = {q_resistance:.1f} kN/m^2")

    print(f"\n--- Final Answer ---")
    print(f"The required ULS design bearing resistance is {required_resistance:.1f} kN/m^2.")

# Run the calculation and print the results
solve_footing_resistance()