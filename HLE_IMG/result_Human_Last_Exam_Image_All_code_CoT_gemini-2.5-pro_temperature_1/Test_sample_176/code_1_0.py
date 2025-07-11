import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation
    based on Eurocode 7 principles for drained conditions.
    """
    # --- Step 1: Define Given Parameters ---
    print("--- Step 1: Given Parameters ---")
    # Loads
    G_k = 500.0  # kN (Permanent vertical load)
    Q_vk = 1500.0  # kN (Variable vertical load)
    Q_hk = 120.0  # kN (Variable horizontal load)

    # Factors (Assuming Eurocode 7, Design Approach 1, Combination 1)
    zeta = 0.7  # Combination factor for independent variable loads
    gamma_G = 1.35  # Partial factor for permanent actions
    gamma_Q = 1.5   # Partial factor for variable actions
    gamma_c_prime = 1.0  # Partial factor for cohesion
    gamma_phi_prime = 1.0 # Partial factor for friction angle
    gamma_gamma = 1.0    # Partial factor for soil unit weight

    # Foundation and Soil Properties
    B = 2.39  # m (Foundation width)
    L = 2.39  # m (Foundation length)
    z = 1.0   # m (Foundation depth)
    h_load = 1.0 # m (Lever arm for horizontal load)
    c_k_prime = 10.0  # kPa (Characteristic effective cohesion)
    phi_k_prime_deg = 28.0  # degrees (Characteristic effective friction angle)
    gamma_s = 20.0  # kN/m^3 (Soil unit weight)

    print(f"Permanent vertical load, G_k = {G_k} kN")
    print(f"Variable vertical load, Q_vk = {Q_vk} kN")
    print(f"Variable horizontal load, Q_hk = {Q_hk} kN")
    print(f"Foundation width, B = L = {B} m")
    print(f"Foundation depth, z = {z} m")
    print(f"Characteristic cohesion, c'_k = {c_k_prime} kPa")
    print(f"Characteristic friction angle, phi'_k = {phi_k_prime_deg} degrees")
    print(f"Soil unit weight, gamma_s = {gamma_s} kN/m^3\n")

    # --- Step 2: Calculate Design Values ---
    print("--- Step 2: Design Values ---")
    # Design Soil Parameters
    c_d_prime = c_k_prime / gamma_c_prime
    phi_d_prime_deg = math.degrees(math.atan(math.tan(math.radians(phi_k_prime_deg)) / gamma_phi_prime))
    phi_d_prime_rad = math.radians(phi_d_prime_deg)
    gamma_d = gamma_s / gamma_gamma
    print(f"Design cohesion, c'_d = {c_d_prime:.2f} kPa")
    print(f"Design friction angle, phi'_d = {phi_d_prime_deg:.2f} degrees")
    print(f"Design soil unit weight, gamma_d = {gamma_d:.2f} kN/m^3")

    # Design Loads (Load Combination 1)
    V_d = gamma_G * G_k + gamma_Q * Q_vk
    H_d = gamma_Q * zeta * Q_hk
    print(f"Design vertical load, V_d = {gamma_G} * {G_k} + {gamma_Q} * {Q_vk} = {V_d:.2f} kN")
    print(f"Design horizontal load, H_d = {gamma_Q} * {zeta} * {Q_hk} = {H_d:.2f} kN\n")

    # --- Step 3: Calculate Effective Foundation Dimensions ---
    print("--- Step 3: Effective Foundation Dimensions ---")
    M_d = H_d * h_load
    e_B = M_d / V_d
    B_prime = B - 2 * e_B
    L_prime = L
    A_prime = B_prime * L_prime
    print(f"Moment at base, M_d = {H_d:.2f} kN * {h_load:.2f} m = {M_d:.2f} kNm")
    print(f"Eccentricity, e_B = {M_d:.2f} / {V_d:.2f} = {e_B:.4f} m")
    print(f"Effective width, B' = {B} - 2 * {e_B:.4f} = {B_prime:.4f} m")
    print(f"Effective length, L' = {L:.4f} m")
    print(f"Effective area, A' = {B_prime:.4f} m * {L_prime:.4f} m = {A_prime:.4f} m^2\n")

    # --- Step 4: Calculate Bearing Capacity Factors ---
    print(f"--- Step 4: Bearing Capacity Factors (for phi'_d = {phi_d_prime_deg:.2f} deg) ---")
    N_q = math.exp(math.pi * math.tan(phi_d_prime_rad)) * (math.tan(math.radians(45) + phi_d_prime_rad / 2))**2
    N_c = (N_q - 1) * (1 / math.tan(phi_d_prime_rad))
    N_gamma = 2 * (N_q - 1) * math.tan(phi_d_prime_rad)
    print(f"N_q = {N_q:.3f}")
    print(f"N_c = {N_c:.3f}")
    print(f"N_gamma = {N_gamma:.3f}\n")

    # --- Step 5: Calculate Shape and Load Inclination Factors ---
    print("--- Step 5: Shape and Load Inclination Factors ---")
    ratio_B_L = B_prime / L_prime
    # Shape factors
    s_q = 1 + (ratio_B_L) * math.sin(phi_d_prime_rad)
    s_gamma = 1 - 0.3 * (ratio_B_L)
    s_c = (s_q * N_q - 1) / (N_q - 1)
    print(f"Shape factor, s_q = {s_q:.4f}")
    print(f"Shape factor, s_gamma = {s_gamma:.4f}")
    print(f"Shape factor, s_c = {s_c:.4f}")
    # Load inclination factors (with simplification from problem statement)
    m = (2 + ratio_B_L) / (1 + ratio_B_L)
    i_q = (1 - H_d / V_d)**m
    i_gamma = (1 - H_d / V_d)**(m + 1)
    i_c = (i_q * N_q - 1) / (N_q - 1)
    print(f"Load inclination factor, i_q = {i_q:.4f}")
    print(f"Load inclination factor, i_gamma = {i_gamma:.4f}")
    print(f"Load inclination factor, i_c = {i_c:.4f}\n")
    
    # --- Step 6: Calculate Overburden and Final Resistance ---
    print("--- Step 6: Final Design Bearing Resistance Calculation ---")
    q_prime = gamma_d * z
    print(f"Overburden pressure, q' = {gamma_d:.2f} kN/m^3 * {z:.2f} m = {q_prime:.2f} kPa\n")
    
    # Base inclination factors (b_c, b_q, b_gamma) are 1.0 as the base is horizontal.
    term1 = c_d_prime * N_c * s_c * i_c
    term2 = q_prime * N_q * s_q * i_q
    term3 = 0.5 * gamma_d * B_prime * N_gamma * s_gamma * i_gamma
    q_ult = term1 + term2 + term3
    R_d = q_ult * A_prime

    print("The formula for design bearing resistance R_d is:")
    print("R_d = A' * (c'_d*N_c*s_c*i_c + q'*N_q*s_q*i_q + 0.5*gamma_d*B'*N_gamma*s_gamma*i_gamma)")
    print("\nPlugging in the calculated values:")
    print(f"R_d = {A_prime:.3f} * ({c_d_prime:.1f} * {N_c:.3f} * {s_c:.3f} * {i_c:.3f} + {q_prime:.1f} * {N_q:.3f} * {s_q:.3f} * {i_q:.3f} + 0.5 * {gamma_d:.1f} * {B_prime:.3f} * {N_gamma:.3f} * {s_gamma:.3f} * {i_gamma:.3f})")
    print(f"R_d = {A_prime:.3f} * ({term1:.2f} + {term2:.2f} + {term3:.2f})")
    print(f"R_d = {A_prime:.3f} * ({q_ult:.2f})")
    print(f"\nFinal Design Resistance, R_d = {R_d:.1f} kN")
    
    return R_d

if __name__ == '__main__':
    design_resistance = calculate_design_resistance()
    # Present the final answer in the requested format
    final_answer = round(design_resistance)
    print(f"\n<<<The final design resistance is approximately {final_answer} kN.>>>")
