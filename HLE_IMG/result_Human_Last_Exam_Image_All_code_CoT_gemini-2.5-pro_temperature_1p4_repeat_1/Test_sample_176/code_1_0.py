import math

def calculate_foundation_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation
    according to Eurocode 7, Load Combination 1 (Drained Conditions).
    """

    # --- 1. Given Parameters ---
    G_k = 500.0  # kN (Permanent vertical load)
    Q_vk = 1500.0  # kN (Variable vertical load)
    Q_hk = 120.0  # kN (Variable horizontal load)
    zeta = 0.7  # Combination factor for independent variable loads
    gamma_c = 24.0  # kN/m^3 (Concrete unit weight)
    Df = 1.0  # m (Foundation depth)
    c_k_prime = 10.0  # kPa (Characteristic effective cohesion)
    phi_k_prime_deg = 28.0  # degrees (Characteristic effective angle of internal friction)
    gamma_soil = 20.0  # kN/m^3 (Soil unit weight)
    B = 2.39  # m (Foundation width)
    L = 2.39  # m (Foundation length)

    # --- 2. Partial Factors (Eurocode 7, DA1/1) ---
    gamma_G = 1.35  # Partial factor for permanent actions
    gamma_Q = 1.50  # Partial factor for variable actions
    gamma_c_prime = 1.0 # Partial factor for cohesion
    gamma_phi_prime = 1.0 # Partial factor for friction angle
    gamma_R_v = 1.0  # Partial factor for bearing resistance

    # --- 3. Calculate Design Actions & Foundation Weight ---
    V_d_structure = gamma_G * G_k + gamma_Q * Q_vk
    H_d = gamma_Q * zeta * Q_hk
    W_f = B * L * Df * gamma_c
    V_d_total = V_d_structure + gamma_G * W_f

    # --- 4. Calculate Design Soil Parameters ---
    phi_k_prime_rad = math.radians(phi_k_prime_deg)
    # Since partial factors are 1.0, design values = characteristic values
    c_d_prime = c_k_prime / gamma_c_prime
    phi_d_prime_deg = math.degrees(math.atan(math.tan(phi_k_prime_rad) / gamma_phi_prime))
    phi_d_prime_rad = math.radians(phi_d_prime_deg)

    # --- 5. Calculate Effective Foundation Dimensions ---
    M_d = H_d * Df
    e = M_d / V_d_total
    B_prime = B - 2 * e
    L_prime = L
    A_prime = B_prime * L_prime

    # --- 6. Calculate Bearing Resistance Components ---
    # a) Overburden pressure
    q = gamma_soil * Df

    # b) Bearing capacity factors
    tan_phi_d = math.tan(phi_d_prime_rad)
    N_q = math.exp(math.pi * tan_phi_d) * (math.tan(math.radians(45) + phi_d_prime_rad / 2) ** 2)
    N_c = (N_q - 1) / tan_phi_d
    N_gamma = 2 * (N_q - 1) * tan_phi_d

    # c) Shape factors
    s_q = 1 + (B_prime / L_prime) * math.sin(phi_d_prime_rad)
    s_gamma = 1 - 0.3 * (B_prime / L_prime)
    s_c = (s_q * N_q - 1) / (N_q - 1)

    # d) Load inclination factors (omitting adhesion term as requested)
    m = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
    i_q = (1 - H_d / V_d_total)**m
    i_gamma = (1 - H_d / V_d_total)**(m + 1)
    i_c = i_q

    # --- 7. Calculate Design Bearing Resistance ---
    term_c = c_d_prime * N_c * s_c * i_c
    term_q = q * N_q * s_q * i_q
    term_gamma = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * i_gamma
    
    R_k = A_prime * (term_c + term_q + term_gamma)
    R_d = R_k / gamma_R_v

    # --- 8. Print Results ---
    print("--- Calculation Steps ---")
    print(f"Design vertical load V_d,total = ({gamma_G} * {G_k} + {gamma_Q} * {Q_vk}) + {gamma_G} * ({B}*{L}*{Df}*{gamma_c}) = {V_d_total:.2f} kN")
    print(f"Design horizontal load H_d = {gamma_Q} * {zeta} * {Q_hk} = {H_d:.2f} kN")
    print(f"Eccentricity e = ({H_d:.2f} * {Df}) / {V_d_total:.2f} = {e:.4f} m")
    print(f"Effective width B' = {B} - 2 * {e:.4f} = {B_prime:.3f} m")
    print(f"Effective area A' = {B_prime:.3f} * {L_prime} = {A_prime:.3f} m²\n")

    print("--- Final Design Resistance Equation ---")
    print("R_d = A' * (c'_d*N_c*s_c*i_c + q*N_q*s_q*i_q + 0.5*γ*B'*N_γ*s_γ*i_γ) / γ_R;v\n")

    print("--- Component Values ---")
    print(f"Bearing Capacity Factors: N_c={N_c:.2f}, N_q={N_q:.2f}, N_γ={N_gamma:.2f}")
    print(f"Shape Factors: s_c={s_c:.3f}, s_q={s_q:.3f}, s_γ={s_gamma:.3f}")
    print(f"Inclination Factors: i_c={i_c:.3f}, i_q={i_q:.3f}, i_γ={i_gamma:.3f}\n")

    print("--- Assembling the Final Equation ---")
    print("R_d =")
    print(f"  {A_prime:.3f} * (")
    print(f"    ({c_d_prime:.1f} * {N_c:.2f} * {s_c:.3f} * {i_c:.3f}) +   // Cohesion term = {term_c:.2f} kPa")
    print(f"    ({q:.1f} * {N_q:.2f} * {s_q:.3f} * {i_q:.3f}) +   // Surcharge term = {term_q:.2f} kPa")
    print(f"    (0.5 * {gamma_soil:.1f} * {B_prime:.3f} * {N_gamma:.2f} * {s_gamma:.3f} * {i_gamma:.3f})   // Self-weight term = {term_gamma:.2f} kPa")
    print(f"  ) / {gamma_R_v}")

    print(f"\nR_d = {A_prime:.3f} * ({term_c:.2f} + {term_q:.2f} + {term_gamma:.2f}) / {gamma_R_v}")
    print(f"R_d = {A_prime:.3f} * {term_c + term_q + term_gamma:.2f} / {gamma_R_v}")
    print(f"\nFinal Design Resistance R_d = {R_d:.1f} kN")
    
    return R_d

# Execute the calculation and print the final answer in the required format
final_answer = calculate_foundation_resistance()
print(f"\n<<<{final_answer:.1f}>>>")
