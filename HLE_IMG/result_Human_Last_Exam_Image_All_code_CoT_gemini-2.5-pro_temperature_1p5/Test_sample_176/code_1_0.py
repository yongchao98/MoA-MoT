import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation
    according to Eurocode 7 (DA1/1) specifications.
    """
    # Step 1: Define Inputs and Partial Factors
    G = 500.0  # Permanent vertical load (kN)
    Qv = 1500.0 # Variable vertical load (kN)
    Qh = 120.0  # Variable horizontal load (kN)
    zeta = 0.7  # Combination factor
    B = 2.39    # Foundation width (m)
    L = 2.39    # Foundation length (m)
    d = 1.0     # Foundation depth (m)
    h_load = 1.0 # Height of horizontal load application above base (m)

    c_k = 10.0  # Characteristic effective cohesion (kPa or kN/m^2)
    phi_k_deg = 28.0 # Characteristic effective friction angle (degrees)
    gamma_soil = 20.0 # Soil unit weight (kN/m^3)
    
    # Partial Factors for Eurocode 7, DA1, Combination 1
    gamma_G = 1.35
    gamma_Q = 1.50
    gamma_c = 1.0
    gamma_phi = 1.0
    gamma_Rv = 1.0 # Partial factor for bearing resistance

    print("--- Input Parameters ---")
    print(f"Permanent Load (G): {G} kN")
    print(f"Variable Vertical Load (Qv): {Qv} kN")
    print(f"Variable Horizontal Load (Qh): {Qh} kN")
    print(f"Foundation Width (B): {B} m, Length (L): {L} m")
    print(f"Foundation Depth (d): {d} m")
    print(f"Soil Cohesion (c'_k): {c_k} kPa, Friction Angle (phi'_k): {phi_k_deg} degrees")
    print("-" * 26 + "\n")

    # Step 2: Calculate Design Actions
    V_d = G * gamma_G + Qv * gamma_Q
    H_d = Qh * gamma_Q * zeta
    
    print("--- Design Actions ---")
    print(f"Design Vertical Load (V_d) = {G} * {gamma_G} + {Qv} * {gamma_Q} = {V_d:.2f} kN")
    print(f"Design Horizontal Load (H_d) = {Qh} * {gamma_Q} * {zeta} = {H_d:.2f} kN")
    print("-" * 22 + "\n")

    # Step 3: Determine Effective Foundation Dimensions
    M_d = H_d * h_load
    e_B = M_d / V_d
    B_prime = B - 2 * e_B
    L_prime = L
    A_prime = B_prime * L_prime

    print("--- Effective Dimensions ---")
    print(f"Eccentricity (e_B) = {M_d:.2f} / {V_d:.2f} = {e_B:.4f} m")
    print(f"Effective Width (B') = {B} - 2 * {e_B:.4f} = {B_prime:.4f} m")
    print(f"Effective Length (L') = {L_prime:.4f} m")
    print(f"Effective Area (A') = {B_prime:.4f} * {L_prime:.4f} = {A_prime:.4f} m^2")
    print("-" * 28 + "\n")

    # Step 4: Calculate Design Soil Properties
    c_d = c_k / gamma_c
    phi_k_rad = math.radians(phi_k_deg)
    phi_d_rad = math.atan(math.tan(phi_k_rad) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)

    print("--- Design Soil Parameters ---")
    print(f"Design Cohesion (c'_d) = {c_k} / {gamma_c} = {c_d:.2f} kPa")
    print(f"Design Friction Angle (phi'_d) = atan(tan({phi_k_deg}) / {gamma_phi}) = {phi_d_deg:.2f} degrees")
    print("-" * 30 + "\n")

    # Step 5: Calculate Bearing Capacity Factors
    N_q = math.exp(math.pi * math.tan(phi_d_rad)) * math.tan(math.radians(45) + phi_d_rad / 2)**2
    N_c = (N_q - 1) / math.tan(phi_d_rad) if math.tan(phi_d_rad) > 1e-6 else (2+math.pi)
    N_gamma = 2 * (N_q - 1) * math.tan(phi_d_rad)

    print("--- Bearing Capacity Factors ---")
    print(f"N_q = {N_q:.2f}")
    print(f"N_c = {N_c:.2f}")
    print(f"N_gamma = {N_gamma:.2f}")
    print("-" * 32 + "\n")

    # Step 6: Calculate Correction Factors
    # Shape factors
    s_q = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
    s_c = (s_q * N_q - 1) / (N_q - 1)
    s_gamma = 1 - 0.3 * (B_prime / L_prime)

    # Depth factors
    d_q = 1 + 2 * math.tan(phi_d_rad) * (1 - math.sin(phi_d_rad))**2 * (d / B_prime)
    d_c = d_q - (1 - d_q) / (N_c * math.tan(phi_d_rad))
    d_gamma = 1.0

    # Inclination factors (simplified as per instruction)
    m = (2 + B_prime / L_prime) / (1 + B_prime / L_prime)
    i_q = (1 - H_d / V_d)**m
    i_gamma = (1 - H_d / V_d)**(m + 1)
    i_c = (i_q * N_q - 1) / (N_q - 1)

    print("--- Correction Factors ---")
    print(f"Shape: s_c={s_c:.3f}, s_q={s_q:.3f}, s_gamma={s_gamma:.3f}")
    print(f"Depth: d_c={d_c:.3f}, d_q={d_q:.3f}, d_gamma={d_gamma:.3f}")
    print(f"Inclination: i_c={i_c:.3f}, i_q={i_q:.3f}, i_gamma={i_gamma:.3f}")
    print("-" * 26 + "\n")
    
    # Step 7: Calculate Design Bearing Resistance
    # Overburden pressure
    q_prime = gamma_soil * d
    
    # Individual terms of the bearing capacity equation
    term_c = c_d * N_c * s_c * d_c * i_c
    term_q = q_prime * N_q * s_q * d_q * i_q
    term_gamma = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * d_gamma * i_gamma

    # Total resistance (characteristic)
    R_k = A_prime * (term_c + term_q + term_gamma)

    # Design resistance
    R_d = R_k / gamma_Rv

    print("--- Design Bearing Resistance Calculation ---")
    print("R_d = A' * (c'_d*N_c*s_c*d_c*i_c + q'*N_q*s_q*d_q*i_q + 0.5*gamma*B'*N_gamma*s_gamma*d_gamma*i_gamma) / gamma_R,v\n")
    print(f"R_d = {A_prime:.3f} * (({c_d:.2f}*{N_c:.2f}*{s_c:.3f}*{d_c:.3f}*{i_c:.3f}) + ({q_prime:.2f}*{N_q:.2f}*{s_q:.3f}*{d_q:.3f}*{i_q:.3f}) + (0.5*{gamma_soil:.2f}*{B_prime:.3f}*{N_gamma:.2f}*{s_gamma:.3f}*{d_gamma:.1f}*{i_gamma:.3f})) / {gamma_Rv:.1f}")
    print(f"R_d = {A_prime:.3f} * ({term_c:.2f} + {term_q:.2f} + {term_gamma:.2f}) / {gamma_Rv:.1f}")
    print(f"R_d = {A_prime:.3f} * ({term_c + term_q + term_gamma:.2f}) / {gamma_Rv:.1f}")
    print(f"R_d = {R_k:.1f} / {gamma_Rv:.1f}")
    print(f"\nFinal Design Resistance (R_d) = {R_d:.1f} kN")
    return R_d

if __name__ == '__main__':
    result = calculate_design_resistance()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<5907.4>>>")
