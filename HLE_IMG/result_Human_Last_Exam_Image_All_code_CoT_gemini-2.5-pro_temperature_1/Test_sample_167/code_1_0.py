import math

def solve_bearing_resistance():
    """
    Calculates the required ULS design bearing resistance for a pad footing.
    """
    # --- 1. Given Parameters ---
    Gk = 1000  # Characteristic permanent vertical load [kN]
    Qk_vert = 1500  # Characteristic variable vertical load [kN]
    Qk_horz = 300  # Characteristic variable horizontal load [kN]
    h_load = 2.0  # Height of horizontal load application above footing base [m]
    Df = 0.75  # Depth of footing base below ground level [m]
    phi_k_deg = 35  # Characteristic angle of internal friction [degrees]
    c_prime_k = 0  # Characteristic effective cohesion [kPa]
    gamma_soil = 20  # Unit weight of soil [kN/m^3]

    # Partial factors for EC7 DA1-1
    gamma_G = 1.35
    gamma_Q = 1.5
    gamma_phi = 1.0
    gamma_c = 1.0
    gamma_gamma = 1.0

    # --- 2. Calculate Design Loads and Soil Properties ---
    # Assumption: Gk includes the self-weight of the foundation.
    Vd = gamma_G * Gk + gamma_Q * Qk_vert
    Hd = gamma_Q * Qk_horz
    Md = Hd * h_load
    e_B = Md / Vd

    # Design soil parameters
    phi_d_rad = math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)
    c_prime_d = c_prime_k / gamma_c
    gamma_prime_d = gamma_soil / gamma_gamma
    
    # Effective overburden pressure at foundation level
    q_prime = gamma_prime_d * Df

    # --- 3. Bearing Capacity Factors (Vesic's formulas, compatible with EC7) ---
    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45 + phi_d_deg / 2)))**2
    N_gamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
    # Nc is not needed as c' = 0

    # --- 4. Iterative Search for Footing Width B ---
    # Find B where applied pressure (sigma_d) equals bearing resistance (q_d)
    B_optimal = 0
    min_abs_diff = float('inf')

    # Iterate through a range of possible B values
    for i in range(100, 5000):
        B = i / 1000.0  # Search B from 0.1m to 5.0m in 1mm steps
        
        # Check for tension: eccentricity must be less than B/6
        if e_B >= B / 6:
            continue
            
        # Effective dimensions
        B_prime = B - 2 * e_B
        L_prime = B  # Since B=L and load is along B-axis
        A_prime = B_prime * L_prime
        
        ratio_B_L = B_prime / L_prime
        
        # Shape factors
        sq = 1 + ratio_B_L * math.sin(phi_d_rad)
        s_gamma = 1 - 0.3 * ratio_B_L
        
        # Inclination factors (for c'=0)
        m = (2 + ratio_B_L) / (1 + ratio_B_L)
        iq = (1 - Hd / Vd)**m
        i_gamma = (1 - Hd / Vd)**(m + 1)
        
        # Design bearing resistance (q_d)
        term_q = q_prime * Nq * sq * iq
        term_gamma = 0.5 * gamma_prime_d * B_prime * N_gamma * s_gamma * i_gamma
        q_d = term_q + term_gamma
        
        # Applied design bearing pressure (sigma_d)
        sigma_d = Vd / A_prime
        
        # Find the B where the difference is minimal
        if abs(q_d - sigma_d) < min_abs_diff:
            min_abs_diff = abs(q_d - sigma_d)
            B_optimal = B

    # --- 5. Calculate Final Resistance with the Optimal B ---
    B = B_optimal
    B_prime = B - 2 * e_B
    L_prime = B
    ratio_B_L = B_prime / L_prime
    sq = 1 + ratio_B_L * math.sin(phi_d_rad)
    s_gamma = 1 - 0.3 * ratio_B_L
    m = (2 + ratio_B_L) / (1 + ratio_B_L)
    iq = (1 - Hd / Vd)**m
    i_gamma = (1 - Hd / Vd)**(m + 1)
    
    term_q = q_prime * Nq * sq * iq
    term_gamma = 0.5 * gamma_prime_d * B_prime * N_gamma * s_gamma * i_gamma
    required_resistance = term_q + term_gamma

    # --- 6. Print Results ---
    print("--- ULS Design Calculation ---")
    print(f"Design Vertical Load (Vd): {Vd:.1f} kN")
    print(f"Design Horizontal Load (Hd): {Hd:.1f} kN")
    print(f"Load Eccentricity (e): {e_B:.3f} m")
    print(f"Optimal Footing Width (B) found: {B:.2f} m")
    print(f"Effective Width (B'): {B_prime:.2f} m\n")

    print("--- Required Design Bearing Resistance Calculation ---")
    print("Formula: q_d = q' * Nq * sq * iq + 0.5 * γ' * B' * Nγ * sγ * iγ\n")
    print("Component values:")
    print(f"q' = {q_prime:.1f} kPa")
    print(f"Nq = {Nq:.2f}")
    print(f"Nγ = {N_gamma:.2f}")
    print(f"sq = {sq:.3f}")
    print(f"sγ = {s_gamma:.3f}")
    print(f"iq = {iq:.3f}")
    print(f"iγ = {i_gamma:.3f}\n")

    print("Final Calculation:")
    print(f"q_d = {q_prime:.1f} * {Nq:.2f} * {sq:.3f} * {iq:.3f} + "
          f"0.5 * {gamma_prime_d:.1f} * {B_prime:.2f} * {N_gamma:.2f} * {s_gamma:.3f} * {i_gamma:.3f}")
    print(f"q_d = {term_q:.1f} + {term_gamma:.1f}")
    print(f"Required ULS design bearing resistance = {required_resistance:.1f} kN/m^2")

solve_bearing_resistance()