import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation
    according to Eurocode 7 principles.
    """

    # --- 1. Define Given Parameters ---
    G_k = 500.0  # Permanent vertical load [kN]
    Q_vk = 1500.0 # Variable vertical load [kN]
    Q_hk = 120.0  # Variable horizontal load [kN]
    zeta = 0.7    # Combination factor for independent variable loads
    B = 2.39      # Foundation width [m]
    L = 2.39      # Foundation length [m]
    D_f = 1.0     # Foundation depth [m]
    c_k = 10.0    # Characteristic effective cohesion [kPa or kN/m^2]
    phi_k_deg = 28.0 # Characteristic effective friction angle [degrees]
    gamma_soil = 20.0 # Soil unit weight [kN/m^3]
    gamma_conc = 24.0 # Concrete unit weight [kN/m^3]

    # Partial factors for Load Combination 1 (DA1-1)
    gamma_G = 1.35
    gamma_Q = 1.5
    gamma_Rv = 1.0 # Partial resistance factor for DA1-1

    # --- 2. Assume Foundation Thickness ---
    # The thickness 'h' is not given, so a reasonable value is assumed.
    h = 0.5  # Assumed foundation thickness [m]
    print(f"Step 1: Input Parameters & Assumptions")
    print(f"Permanent Load G_k = {G_k} kN")
    print(f"Variable Vertical Load Q_vk = {Q_vk} kN")
    print(f"Variable Horizontal Load Q_hk = {Q_hk} kN")
    print(f"Foundation Width B = {B} m")
    print(f"Foundation Depth D_f = {D_f} m")
    print(f"Assumed Foundation Thickness h = {h} m")
    print(f"Soil properties: c'_k = {c_k} kPa, φ'_k = {phi_k_deg}°, γ_soil = {gamma_soil} kN/m³")
    print("-" * 30)

    # --- 3. Calculate Design Actions ---
    # Characteristic self-weight of the foundation
    G_found = B * L * h * gamma_conc
    # Design vertical load
    V_d = gamma_G * (G_k + G_found) + gamma_Q * Q_vk
    # Design horizontal load
    H_d = gamma_Q * zeta * Q_hk
    # Design moment at the base (H_d applied 1.0m above base)
    M_d = H_d * D_f
    print(f"Step 2: Design Actions (Loads)")
    print(f"Foundation self-weight G_found = {G_found:.2f} kN")
    print(f"Design Vertical Load V_d = 1.35 * ({G_k} + {G_found:.2f}) + 1.5 * {Q_vk} = {V_d:.2f} kN")
    print(f"Design Horizontal Load H_d = 1.5 * {zeta} * {Q_hk} = {H_d:.2f} kN")
    print(f"Design Moment M_d = {H_d:.2f} * {D_f} = {M_d:.2f} kNm")
    print("-" * 30)
    
    # --- 4. Determine Effective Foundation Dimensions ---
    e = M_d / V_d
    B_prime = B - 2 * e
    L_prime = L  # No eccentricity in L direction
    A_prime = B_prime * L_prime
    print(f"Step 3: Effective Foundation Dimensions")
    print(f"Eccentricity e = {M_d:.2f} / {V_d:.2f} = {e:.4f} m")
    print(f"Effective Width B' = {B} - 2 * {e:.4f} = {B_prime:.4f} m")
    print(f"Effective Length L' = {L_prime:.4f} m")
    print(f"Effective Area A' = {B_prime:.4f} * {L_prime:.4f} = {A_prime:.4f} m²")
    print("-" * 30)
    
    # --- 5. Calculate Bearing Resistance Factors ---
    phi_k_rad = math.radians(phi_k_deg)

    # Bearing Capacity Factors (EC7 Annex D)
    Nq = math.exp(math.pi * math.tan(phi_k_rad)) * math.tan(math.radians(45 + phi_k_deg / 2))**2
    Nc = (Nq - 1) / math.tan(phi_k_rad)
    N_gamma = 2 * (Nq - 1) * math.tan(phi_k_rad)

    # Shape Factors
    ratio_BL = B_prime / L_prime
    sq = 1 + ratio_BL * math.sin(phi_k_rad)
    s_gamma = 1 - 0.3 * ratio_BL
    sc = (sq * Nq - 1) / (Nq - 1)
    
    # Inclination Factors (with simplification from prompt)
    m = (2 + ratio_BL) / (1 + ratio_BL)
    ratio_HV = H_d / V_d
    # Per prompt, omit A'c'cot(phi) term
    iq = (1 - ratio_HV)**m
    i_gamma = (1 - ratio_HV)**(m + 1)
    # i_c is dependent on i_q
    ic = iq - (1 - iq) / (Nc * math.tan(phi_k_rad))

    print(f"Step 4: Bearing Capacity, Shape, and Inclination Factors")
    print(f"Bearing Capacity Factors: N_c = {Nc:.2f}, N_q = {Nq:.2f}, N_γ = {N_gamma:.2f}")
    print(f"Shape Factors: s_c = {sc:.3f}, s_q = {sq:.3f}, s_γ = {s_gamma:.3f}")
    print(f"Inclination Factors: i_c = {ic:.3f}, i_q = {iq:.3f}, i_γ = {i_gamma:.3f}")
    print("-" * 30)
    
    # --- 6. Calculate Design Bearing Resistance ---
    # Design soil parameters are same as characteristic for DA1-1
    c_d = c_k
    gamma_d = gamma_soil

    # Overburden pressure at foundation base
    q_prime = gamma_d * D_f
    
    # Bearing capacity equation terms
    term_c = c_d * Nc * sc * ic
    term_q = q_prime * Nq * sq * iq
    term_gamma = 0.5 * gamma_d * B_prime * N_gamma * s_gamma * i_gamma
    
    # Characteristic resistance (per unit area)
    q_ult = term_c + term_q + term_gamma
    
    # Characteristic resistance (total)
    R_k = q_ult * A_prime
    
    # Design Resistance
    R_d = R_k / gamma_Rv

    # --- 7. Output Results ---
    print("Step 5: Calculate Final Design Resistance R_d")
    print("The final equation is R_d = (Term_c + Term_q + Term_γ) * A' / γ_R,v")
    print("Where:")
    print(f"  Term_c = c'_d * N_c * s_c * i_c = {c_d:.1f} * {Nc:.2f} * {sc:.3f} * {ic:.3f} = {term_c:.2f} kPa")
    print(f"  Term_q = q' * N_q * s_q * i_q = {q_prime:.1f} * {Nq:.2f} * {sq:.3f} * {iq:.3f} = {term_q:.2f} kPa")
    print(f"  Term_γ = 0.5 * γ' * B' * N_γ * s_γ * i_γ = 0.5 * {gamma_d:.1f} * {B_prime:.3f} * {N_gamma:.2f} * {s_gamma:.3f} * {i_gamma:.3f} = {term_gamma:.2f} kPa")
    print(f"Characteristic bearing pressure q_ult = {term_c:.2f} + {term_q:.2f} + {term_gamma:.2f} = {q_ult:.2f} kPa")
    print("-" * 30)
    print("Final Calculation:")
    print(f"Design Resistance R_d = ({q_ult:.2f} kPa * {A_prime:.3f} m²) / {gamma_Rv}")
    print(f"Design Resistance R_d = {R_d:.2f} kN")

calculate_design_resistance()
<<<5355.26>>>