import math

def solve_bearing_resistance():
    """
    Calculates the required ULS design bearing resistance for a pad footing
    by iteratively finding the footing size where design action equals design resistance.
    """
    # 1. Define constants and parameters from the problem
    # Loads
    Gk = 1000.0  # kN (Permanent vertical load)
    Qk_v = 1500.0 # kN (Variable vertical load)
    Qk_h = 300.0  # kN (Variable horizontal load)

    # Geometry
    h_f = 0.75  # m (Footing thickness, assumed equal to depth D)
    D = 0.75    # m (Depth of footing base)
    h_col = 2.0 # m (Height of column for moment arm)

    # Material Properties
    gamma_c = 24.0  # kN/m^3 (Unit weight of concrete)
    gamma_s = 20.0  # kN/m^3 (Unit weight of soil)
    phi_k_deg = 35.0 # degrees (Characteristic angle of shearing resistance)
    c_k = 0.0       # kPa (Characteristic effective cohesion)

    # Eurocode 7 Partial Factors (Design Approach 1, Combination 2)
    gamma_G = 1.0   # Permanent actions
    gamma_Q = 1.3   # Variable actions
    gamma_phi = 1.25 # Angle of shearing resistance
    gamma_c_prime = 1.25 # Effective cohesion (not used as c'=0)
    gamma_R_v = 1.4 # Bearing resistance

    # 2. Setup for Iterative Search
    B = 4.0  # m (Start with a safe guess for B)
    step = 0.001 # m (Iteration step size for precision)
    
    # Pre-calculate design soil parameters
    phi_k_rad = math.radians(phi_k_deg)
    phi_d_rad = math.atan(math.tan(phi_k_rad) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)
    c_d = c_k / gamma_c_prime

    # Pre-calculate Bearing Capacity Factors for design friction angle phi_d
    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45 + phi_d_deg / 2)))**2
    Nc = (Nq - 1) / math.tan(phi_d_rad) if math.tan(phi_d_rad) > 1e-6 else 0
    N_gamma = 2 * (Nq - 1) * math.tan(phi_d_rad)

    # Surcharge pressure at foundation level
    q_prime = gamma_s * D
    
    last_safe_values = {}

    # Loop to find B where q_Ed is just less than or equal to q_Rd
    while B > 0:
        # Calculate ULS Design Actions for the current footing size B
        W_f = B * B * h_f * gamma_c
        Vd = gamma_G * (Gk + W_f) + gamma_Q * Qk_v
        Hd = gamma_Q * Qk_h
        Md = Hd * (h_col + h_f)

        # Calculate eccentricity and effective dimensions
        if Vd <= 0: break
        e = Md / Vd
        if e >= B / 2: # No contact area, invalid size
            B -= step
            continue

        B_prime = B - 2 * e
        L_prime = B 
        A_prime = B_prime * L_prime

        if A_prime <= 0:
            B -= step
            continue

        # Calculate Design Bearing Pressure (Action), q_Ed
        q_Ed = Vd / A_prime

        # Calculate Design Bearing Resistance (Capacity), q_Rd
        # Shape factors
        sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
        s_gamma = 1 - 0.3 * (B_prime / L_prime)
        
        # Depth factors
        dq = 1 + 2 * math.tan(phi_d_rad) * (1 - math.sin(phi_d_rad))**2 * (D/B)
        d_gamma = 1
        
        # Inclination factors
        m_exp = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
        # Since c_d is 0, the denominator simplifies
        iq = (1 - Hd / Vd)**m_exp
        i_gamma = (1 - Hd / Vd)**(m_exp + 1)
        
        # Ultimate bearing capacity q_ult (c' term is zero)
        term_q = q_prime * Nq * sq * dq * iq
        term_gamma = 0.5 * gamma_s * B_prime * N_gamma * s_gamma * d_gamma * i_gamma
        q_ult = term_q + term_gamma

        # Design bearing resistance q_Rd
        q_Rd = q_ult / gamma_R_v

        if q_Ed <= q_Rd:
            # Store the values of this valid iteration
            last_safe_values = {'B': B, 'Vd': Vd, 'Hd': Hd, 'Md': Md, 'e': e, 'B_prime': B_prime, 'A_prime': A_prime, 'q_Ed': q_Ed, 'q_Rd': q_Rd}
            B -= step
        else:
            # q_Ed has exceeded q_Rd, so the previous step was the optimal solution
            break

    # 3. Print the final results
    required_resistance = last_safe_values.get('q_Ed', 0)
    
    print("This problem requires finding the footing size 'B' where the applied design pressure (q_Ed) equals the soil's design bearing resistance (q_Rd).")
    print("\n--- Final Design Parameters at Limit State ---")
    print(f"Optimal Footing Width B = {last_safe_values.get('B', 'N/A'):.2f} m")
    print(f"Design Vertical Load Vd = {last_safe_values.get('Vd', 'N/A'):.2f} kN")
    print(f"Design Horizontal Load Hd = {gamma_Q:.1f} * {Qk_h:.0f} = {last_safe_values.get('Hd', 'N/A'):.2f} kN")
    print(f"Design Moment Md = {last_safe_values.get('Hd', 'N/A'):.2f} kN * {h_col + h_f:.2f} m = {last_safe_values.get('Md', 'N/A'):.2f} kNm")
    print(f"Eccentricity e = {last_safe_values.get('e', 'N/A'):.4f} m")
    print(f"Effective Width B' = {last_safe_values.get('B', 'N/A'):.2f} m - 2 * {last_safe_values.get('e', 'N/A'):.4f} m = {last_safe_values.get('B_prime', 'N/A'):.4f} m")
    print(f"Effective Area A' = {last_safe_values.get('B_prime', 'N/A'):.4f} m * {last_safe_values.get('B', 'N/A'):.2f} m = {last_safe_values.get('A_prime', 'N/A'):.4f} m^2")
    
    print("\n--- Final Calculation ---")
    print("The required ULS design bearing resistance is the design pressure at this limit state.")
    print(f"q_required = Vd / A' = {last_safe_values.get('Vd', 'N/A'):.2f} kN / {last_safe_values.get('A_prime', 'N/A'):.4f} m^2")
    print(f"Required ULS Design Bearing Resistance = {required_resistance:.1f} kN/m^2")

solve_bearing_resistance()