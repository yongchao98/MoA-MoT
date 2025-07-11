import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation
    according to Eurocode 7, Design Approach 1, Combination 1.
    """
    # --- Given Parameters ---
    G = 500.0  # Permanent vertical load, kN
    Q_v = 1500.0  # Variable vertical load, kN
    Q_h = 120.0  # Variable horizontal load, kN
    zeta = 0.7  # Combination factor
    D = 1.0  # Foundation depth, m
    c_k = 10.0  # Characteristic effective cohesion, kPa (kN/m^2)
    phi_k_deg = 28.0  # Characteristic effective friction angle, degrees
    gamma_s = 20.0  # Soil unit weight, kN/m^3
    B = 2.39  # Foundation width, m
    L = 2.39  # Foundation length, m
    h_lever_arm = 1.0  # Lever arm for horizontal load (assumed foundation thickness), m

    # --- Eurocode 7, DA1-1 Partial Factors ---
    gamma_G = 1.35
    gamma_Q = 1.50

    # 1. Calculate Design Loads
    V_d = gamma_G * G + gamma_Q * Q_v
    H_d = gamma_Q * zeta * Q_h

    # 2. Calculate Eccentricity and Effective Dimensions
    M_d = H_d * h_lever_arm
    e_B = M_d / V_d
    B_prime = B - 2 * e_B
    L_prime = L  # No eccentricity in L direction
    A_prime = B_prime * L_prime

    # 3. Calculate Overburden Pressure
    q_prime = gamma_s * D

    # 4. Calculate Bearing Capacity Factors
    phi_k_rad = math.radians(phi_k_deg)
    tan_phi = math.tan(phi_k_rad)
    
    N_q = math.exp(math.pi * tan_phi) * math.pow(math.tan(math.radians(45) + phi_k_rad / 2), 2)
    # Check for phi = 0 case
    if phi_k_deg > 0:
        N_c = (N_q - 1) / tan_phi
    else: # Cohesion-only case (undrained)
        N_c = 5.14 
    N_gamma = 2 * (N_q - 1) * tan_phi

    # 5. Calculate Shape Factors
    sin_phi = math.sin(phi_k_rad)
    s_q = 1 + (B_prime / L_prime) * sin_phi
    s_gamma = 1 - 0.3 * (B_prime / L_prime)
    s_c = (s_q * N_q - 1) / (N_q - 1)

    # 6. Calculate Load Inclination Factors
    m = (2 + B_prime / L_prime) / (1 + B_prime / L_prime)
    # Per instruction, omitting the A'*c'*cot(phi) term in the denominator
    i_q = math.pow(1 - H_d / V_d, m)
    i_c = i_q
    i_gamma = math.pow(1 - H_d / V_d, m + 1)

    # 7. Calculate Design Bearing Resistance (R_d)
    # In DA1-1, partial factor for resistance is 1.0, so R_d = R_k
    term_c = c_k * N_c * s_c * i_c
    term_q = q_prime * N_q * s_q * i_q
    term_gamma = 0.5 * gamma_s * B_prime * N_gamma * s_gamma * i_gamma
    
    # Ultimate bearing pressure q_ult in kPa
    q_ult = term_c + term_q + term_gamma
    
    # Design resistance R_d in kN
    R_d = q_ult * A_prime
    
    # --- Print Results ---
    print("--- Calculation of Design Bearing Resistance (R_d) ---")
    
    print("\n1. Final Equation for Design Resistance R_d:")
    print("R_d = A' * (c'_k * N_c * s_c * i_c + q' * N_q * s_q * i_q + 0.5 * γ' * B' * N_γ * s_γ * i_γ)")
    
    print(f"\nR_d = {A_prime:.3f} * ({c_k} * {N_c:.3f} * {s_c:.3f} * {i_c:.3f} + "
          f"{q_prime} * {N_q:.3f} * {s_q:.3f} * {i_q:.3f} + "
          f"0.5 * {gamma_s} * {B_prime:.3f} * {N_gamma:.3f} * {s_gamma:.3f} * {i_gamma:.3f})")
    
    print(f"\nR_d = {A_prime:.3f} * ({term_c:.2f} [kPa] + {term_q:.2f} [kPa] + {term_gamma:.2f} [kPa])")
    print(f"\nR_d = {A_prime:.3f} [m^2] * {q_ult:.2f} [kPa]")
    
    print("\n--- Final Answer ---")
    print(f"The design resistance under load combination 1 is: {R_d:.1f} kN")

calculate_design_resistance()