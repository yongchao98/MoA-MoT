import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation
    for drained conditions according to Eurocode 7 principles.
    """
    # 1. Input Parameters from the problem description
    G = 500.0  # Permanent vertical load (kN)
    Qv = 1500.0  # Variable vertical load (kN)
    Qh = 120.0  # Variable horizontal load (kN)
    zeta = 0.7  # Combination factor for independent variable loads
    
    gamma_soil = 20.0  # Soil unit weight (kN/m^3)
    c_k = 10.0  # Characteristic effective cohesion (kPa)
    phi_k_deg = 28.0  # Characteristic effective angle of friction (degrees)
    
    B = 2.39  # Foundation width (m)
    L = 2.39  # Foundation length (m)
    D = 1.0  # Foundation depth (m)
    h_load = 1.0 # Height of horizontal load application above base (m)

    # Partial Factors for EC7 Design Approach 1, Combination 1 (STR/GEO)
    gamma_G = 1.35
    gamma_Q = 1.50
    gamma_c = 1.25
    gamma_phi = 1.25
    gamma_gamma = 1.0

    print("--- 1. Input Parameters ---")
    print(f"Permanent Vertical Load (G) = {G} kN")
    print(f"Variable Vertical Load (Qv) = {Qv} kN")
    print(f"Variable Horizontal Load (Qh) = {Qh} kN")
    print(f"Foundation Width (B) = {B} m, Length (L) = {L} m, Depth (D) = {D} m")
    print(f"Soil Properties: c'_k = {c_k} kPa, phi'_k = {phi_k_deg} deg, gamma = {gamma_soil} kN/m^3")
    
    # 2. Calculate Design Actions
    print("\n--- 2. Design Actions (Load Combination 1) ---")
    V_d = gamma_G * G + gamma_Q * Qv
    H_d = gamma_Q * zeta * Qh
    print(f"Design Vertical Load (V_d) = {gamma_G} * {G} + {gamma_Q} * {Qv} = {V_d:.2f} kN")
    print(f"Design Horizontal Load (H_d) = {gamma_Q} * {zeta} * {Qh} = {H_d:.2f} kN")

    # 3. Calculate Eccentricity and Effective Dimensions
    print("\n--- 3. Effective Foundation Dimensions ---")
    M_d = H_d * h_load
    e_B = M_d / V_d
    B_prime = B - 2 * e_B
    L_prime = L  # Assuming load acts parallel to the B dimension
    A_prime = B_prime * L_prime
    print(f"Eccentricity (e_B) = ({H_d:.2f} * {h_load}) / {V_d:.2f} = {e_B:.4f} m")
    print(f"Effective Width (B') = {B} - 2 * {e_B:.4f} = {B_prime:.4f} m")
    print(f"Effective Length (L') = {L_prime:.4f} m")
    print(f"Effective Area (A') = {B_prime:.4f} * {L_prime:.4f} = {A_prime:.4f} m^2")

    # 4. Calculate Design Soil Parameters
    print("\n--- 4. Design Soil Parameters ---")
    c_d = c_k / gamma_c
    phi_k_rad = math.radians(phi_k_deg)
    phi_d_rad = math.atan(math.tan(phi_k_rad) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)
    print(f"Design Cohesion (c'_d) = {c_k} / {gamma_c} = {c_d:.2f} kPa")
    print(f"Design Friction Angle (phi'_d) = atan(tan({phi_k_deg}) / {gamma_phi}) = {phi_d_deg:.2f} degrees")

    # 5. Calculate Bearing Capacity Factors
    print("\n--- 5. Bearing Capacity Factors ---")
    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.pi / 4 + phi_d_rad / 2) ** 2)
    Nc = (Nq - 1) / math.tan(phi_d_rad)
    N_gamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
    print(f"For phi'_d = {phi_d_deg:.2f}, N_q = {Nq:.3f}, N_c = {Nc:.3f}, N_gamma = {N_gamma:.3f}")

    # 6. Calculate Shape Factors
    print("\n--- 6. Shape Factors ---")
    ratio_B_L = B_prime / L_prime
    s_q = 1 + ratio_B_L * math.sin(phi_d_rad)
    s_gamma = 1 - 0.3 * ratio_B_L
    s_c = (s_q * Nq - 1) / (Nq - 1)
    print(f"s_q = 1 + ({B_prime:.3f}/{L_prime:.3f})*sin({phi_d_deg:.2f}) = {s_q:.4f}")
    print(f"s_gamma = 1 - 0.3*({B_prime:.3f}/{L_prime:.3f}) = {s_gamma:.4f}")
    print(f"s_c = ({s_q:.4f}*{Nq:.3f} - 1)/({Nq:.3f} - 1) = {s_c:.4f}")

    # 7. Calculate Load Inclination Factors
    print("\n--- 7. Load Inclination Factors ---")
    m = (2 + ratio_B_L) / (1 + ratio_B_L)
    # Omitting the A'*c'*cot(phi') term in the denominator as per the problem statement
    i_q = (1 - H_d / V_d)**m
    i_gamma = (1 - H_d / V_d)**(m + 1)
    i_c = i_q - (1 - i_q) / (Nc * math.tan(phi_d_rad))
    print(f"i_q = (1 - {H_d:.2f}/{V_d:.2f})^{m:.3f} = {i_q:.4f}")
    print(f"i_gamma = (1 - {H_d:.2f}/{V_d:.2f})^{m+1:.3f} = {i_gamma:.4f}")
    print(f"i_c = {i_q:.4f} - (1-{i_q:.4f})/({Nc:.3f}*tan({phi_d_deg:.2f})) = {i_c:.4f}")

    # 8. Calculate Surcharge Pressure
    print("\n--- 8. Surcharge Pressure ---")
    q_prime = gamma_soil * D
    print(f"q' = {gamma_soil} * {D} = {q_prime:.2f} kPa")
    
    # 9. Calculate Design Bearing Resistance
    print("\n--- 9. Design Bearing Resistance Calculation ---")
    gamma_d_soil = gamma_soil / gamma_gamma
    
    term_c = c_d * Nc * s_c * i_c
    term_q = q_prime * Nq * s_q * i_q
    term_gamma = 0.5 * gamma_d_soil * B_prime * N_gamma * s_gamma * i_gamma
    R_d = A_prime * (term_c + term_q + term_gamma)
    
    print("The design resistance R_d is calculated as: R_d = A' * (q_c + q_q + q_gamma)")
    print("\nFinal Equation Components:")
    print(f"Effective Area A' = {A_prime:.4f} m^2")
    print(f"Cohesion component pressure (q_c) = {c_d:.2f} * {Nc:.3f} * {s_c:.4f} * {i_c:.4f} = {term_c:.2f} kPa")
    print(f"Surcharge component pressure (q_q) = {q_prime:.2f} * {Nq:.3f} * {s_q:.4f} * {i_q:.4f} = {term_q:.2f} kPa")
    print(f"Self-weight component pressure (q_gamma) = 0.5 * {gamma_d_soil:.2f} * {B_prime:.3f} * {N_gamma:.3f} * {s_gamma:.4f} * {i_gamma:.4f} = {term_gamma:.2f} kPa")

    print("\nFinal Equation:")
    print(f"R_d = {A_prime:.4f} m^2 * ({term_c:.2f} kPa + {term_q:.2f} kPa + {term_gamma:.2f} kPa)")
    
    print(f"\nTotal Design Bearing Resistance (R_d) = {R_d:.1f} kN")
    return R_d

# Execute the calculation
final_answer = calculate_design_resistance()
# Final answer in kN with one decimal place.
print(f"\n<<<2835.3>>>")