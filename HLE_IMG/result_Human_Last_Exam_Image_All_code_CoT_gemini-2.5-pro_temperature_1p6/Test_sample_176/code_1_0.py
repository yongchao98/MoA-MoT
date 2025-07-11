import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation.
    """
    # 1. Define Inputs
    G_k = 500.0  # Permanent vertical load, kN
    Q_v_k = 1500.0 # Variable vertical load, kN
    Q_h_k = 120.0  # Variable horizontal load, kN
    zeta = 0.7   # Combination factor for independent variable loads
    
    B = 2.39  # Foundation width, m
    L = 2.39  # Foundation length, m
    d = 1.0   # Foundation depth, m
    # Assuming the foundation thickness 'h' is equal to the depth 'd'
    # as the horizontal load is at the top of the foundation.
    h = 1.0   # Foundation thickness, m
    
    gamma_conc = 24.0  # Unit weight of concrete, kN/m^3
    gamma_soil = 20.0  # Unit weight of soil, kN/m^3
    c_k_prime = 10.0   # Characteristic effective cohesion, kPa (kN/m^2)
    phi_k_prime = 28.0 # Characteristic effective angle of shearing resistance, degrees
    
    # Partial factors for actions (EC7 DA1-1)
    gamma_G = 1.35
    gamma_Q = 1.5
    
    # 2. Determine Design Actions (Loads)
    # Foundation self-weight (permanent action)
    W_f = B * L * h * gamma_conc
    G_k_total = G_k + W_f
    
    # Design vertical load
    V_d = gamma_G * G_k_total + gamma_Q * Q_v_k
    # Design horizontal load (assuming vertical is leading, horizontal is accompanying)
    H_d = gamma_Q * zeta * Q_h_k
    
    # 3. Determine Design Soil Parameters (EC7 DA1-1, Set M1)
    gamma_phi_prime = 1.0
    gamma_c_prime = 1.0
    gamma_gamma = 1.0
    
    c_d_prime = c_k_prime / gamma_c_prime
    phi_d_prime_rad = math.atan(math.tan(math.radians(phi_k_prime)) / gamma_phi_prime)
    phi_d_prime_deg = math.degrees(phi_d_prime_rad)
    gamma_d_soil = gamma_soil / gamma_gamma
    
    # Overburden pressure at foundation base
    q = gamma_d_soil * d
    
    # 4. Calculate Effective Foundation Dimensions
    # Moment at the base
    M_d = H_d * h
    # Eccentricity
    e = M_d / V_d
    # Effective width
    B_prime = B - 2 * e
    L_prime = L
    A_prime = B_prime * L_prime
    
    # 5. Calculate Bearing Capacity Factors
    tan_phi = math.tan(phi_d_prime_rad)
    sin_phi = math.sin(phi_d_prime_rad)
    
    N_q = math.exp(math.pi * tan_phi) * (math.tan(math.radians(45) + phi_d_prime_rad / 2))**2
    N_c = (N_q - 1) / tan_phi
    N_gamma = 2 * (N_q - 1) * tan_phi # Vesic formula
    
    # 6. Calculate Correction Factors
    # Shape factors
    ratio_B_L = B_prime / L_prime
    s_q = 1 + ratio_B_L * sin_phi
    s_gamma = 1 - 0.3 * ratio_B_L
    s_c = (s_q * N_q - 1) / (N_q - 1)
    
    # Depth factors
    k = d / B if d / B <= 1 else 1 # Hansen's k=d/B for d/B <=1
    d_q = 1 + 2 * tan_phi * (1 - sin_phi)**2 * k
    d_gamma = 1.0
    d_c = d_q - (1 - d_q) / (N_c * tan_phi)

    # Inclination factors
    # Assuming horizontal load acts parallel to width B
    m = (2 + ratio_B_L) / (1 + ratio_B_L)
    # The problem specifies to omit the A'c'cot(phi) term
    # V_d in the denominator is the total vertical load on the base.
    i_q = (1 - H_d / V_d)**m
    i_gamma = (1 - H_d / V_d)**(m + 1)
    i_c = i_q - (1 - i_q) / (N_c * tan_phi)

    # 7. Calculate Design Bearing Resistance
    term_c = c_d_prime * N_c * s_c * d_c * i_c
    term_q = q * N_q * s_q * d_q * i_q
    term_gamma = 0.5 * gamma_d_soil * B_prime * N_gamma * s_gamma * d_gamma * i_gamma
    
    R_d = A_prime * (term_c + term_q + term_gamma)
    
    # 8. Output the results
    print("--- Intermediate Values ---")
    print(f"Foundation Self-Weight W_f = {W_f:.2f} kN")
    print(f"Total Permanent Load G_k,total = {G_k_total:.2f} kN")
    print(f"Design Vertical Load V_d = {V_d:.2f} kN")
    print(f"Design Horizontal Load H_d = {H_d:.2f} kN")
    print(f"Eccentricity e = {e:.4f} m")
    print(f"Effective Width B' = {B_prime:.4f} m")
    print(f"Effective Area A' = {A_prime:.4f} m^2")
    print(f"Bearing Capacity Factors: N_c = {N_c:.3f}, N_q = {N_q:.3f}, N_gamma = {N_gamma:.3f}")
    print(f"Shape Factors: s_c = {s_c:.3f}, s_q = {s_q:.3f}, s_gamma = {s_gamma:.3f}")
    print(f"Depth Factors: d_c = {d_c:.3f}, d_q = {d_q:.3f}, d_gamma = {d_gamma:.3f}")
    print(f"Inclination Factors: i_c = {i_c:.3f}, i_q = {i_q:.3f}, i_gamma = {i_gamma:.3f}")
    print("\n--- Final Calculation ---")
    print("The design resistance R_d is calculated using the formula:")
    print("R_d = A' * (c'_d*N_c*s_c*d_c*i_c + q*N_q*s_q*d_q*i_q + 0.5*γ'*B'*N_γ*s_γ*d_γ*i_γ)")
    print("\nSubstituting the values:")
    print(f"R_d = {A_prime:.4f} * ({c_d_prime:.2f}*{N_c:.3f}*{s_c:.3f}*{d_c:.3f}*{i_c:.3f} + "
          f"{q:.2f}*{N_q:.3f}*{s_q:.3f}*{d_q:.3f}*{i_q:.3f} + "
          f"0.5*{gamma_d_soil:.2f}*{B_prime:.4f}*{N_gamma:.3f}*{s_gamma:.3f}*{d_gamma:.1f}*{i_gamma:.3f})")
    
    print(f"R_d = {A_prime:.4f} * ({term_c:.2f} + {term_q:.2f} + {term_gamma:.2f}) kN/m^2")
    print(f"R_d = {A_prime:.4f} * ({term_c + term_q + term_gamma:.2f}) kN/m^2")
    print(f"\nThe calculated Design Bearing Resistance is:\n R_d = {R_d:.1f} kN")

    return R_d

if __name__ == "__main__":
    design_resistance = calculate_design_resistance()
    # The final answer in the specified format
    print(f"\n<<<{design_resistance:.1f}>>>")
