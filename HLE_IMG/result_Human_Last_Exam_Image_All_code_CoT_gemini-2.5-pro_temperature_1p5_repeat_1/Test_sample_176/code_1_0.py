import math

def calculate_design_resistance():
    """
    Calculates the design bearing resistance of a square pad foundation.
    """

    # 1. Given Data and Partial Factors
    # Loads
    G_k = 500.0  # kN (Permanent vertical load)
    Q_v_k = 1500.0 # kN (Variable vertical load)
    Q_h_k = 120.0  # kN (Variable horizontal load)

    # Factors
    zeta = 0.7        # Combination factor
    gamma_G = 1.35    # Partial factor for permanent action
    gamma_Q = 1.50    # Partial factor for variable action
    gamma_c = 1.25    # Partial factor for cohesion
    gamma_phi = 1.25  # Partial factor for friction angle

    # Foundation and Soil Properties
    B = 2.39          # m (Foundation width)
    L = 2.39          # m (Foundation length)
    D = 1.0           # m (Foundation depth)
    c_k_prime = 10.0  # kPa (Characteristic effective cohesion)
    phi_k_prime_deg = 28.0 # degrees (Characteristic effective friction angle)
    gamma_soil = 20.0 # kN/m^3 (Soil unit weight)

    print("Step 1: Calculate Design Actions and Soil Parameters")
    # Design Loads
    V_d = gamma_G * G_k + gamma_Q * Q_v_k
    H_d = gamma_Q * (zeta * Q_h_k)
    print(f"Design Vertical Load (Vd) = {gamma_G} * {G_k} + {gamma_Q} * {Q_v_k} = {V_d:.2f} kN")
    print(f"Design Horizontal Load (Hd) = {gamma_Q} * ({zeta} * {Q_h_k}) = {H_d:.2f} kN")

    # Design Soil Parameters
    c_d_prime = c_k_prime / gamma_c
    phi_k_prime_rad = math.radians(phi_k_prime_deg)
    phi_d_prime_rad = math.atan(math.tan(phi_k_prime_rad) / gamma_phi)
    phi_d_prime_deg = math.degrees(phi_d_prime_rad)
    print(f"Design Cohesion (c'd) = {c_k_prime} / {gamma_c} = {c_d_prime:.2f} kPa")
    print(f"Design Friction Angle (φ'd) = arctan(tan({phi_k_prime_deg}°) / {gamma_phi}) = {phi_d_prime_deg:.2f}°")
    print("-" * 30)

    print("Step 2: Calculate Effective Foundation Dimensions")
    # Moment M = Hd * D (assuming horizontal load acts at ground level, h=D)
    M = H_d * D
    e = M / V_d
    B_prime = B - 2 * e
    L_prime = L
    A_prime = B_prime * L_prime
    print(f"Moment (M) = {H_d:.2f} kN * {D:.1f} m = {M:.2f} kNm")
    print(f"Eccentricity (e) = {M:.2f} / {V_d:.2f} = {e:.4f} m")
    print(f"Effective Width (B') = {B:.2f} - 2 * {e:.4f} = {B_prime:.4f} m")
    print(f"Effective Length (L') = {L:.2f} m")
    print(f"Effective Area (A') = {B_prime:.4f} * {L_prime:.2f} = {A_prime:.4f} m²")
    print("-" * 30)

    print("Step 3: Calculate Bearing Capacity Factors (for φ'd = {:.2f}°)".format(phi_d_prime_deg))
    tan_phi_d = math.tan(phi_d_prime_rad)
    N_q = math.exp(math.pi * tan_phi_d) * (math.tan(math.radians(45) + phi_d_prime_rad / 2)**2)
    N_c = (N_q - 1) / tan_phi_d if tan_phi_d > 0 else 5.14
    N_gamma = 2 * (N_q - 1) * tan_phi_d
    print(f"Nq = {N_q:.2f}")
    print(f"Nc = {N_c:.2f}")
    print(f"Nγ = {N_gamma:.2f}")
    print("-" * 30)
    
    print("Step 4: Calculate Shape, Depth, and Inclination Factors")
    # Shape Factors
    B_L_ratio = B_prime / L_prime
    s_q = 1 + B_L_ratio * math.sin(phi_d_prime_rad)
    s_gamma = 1 - 0.3 * B_L_ratio
    s_c = (s_q * N_q - 1) / (N_q - 1)
    print(f"Shape Factors: sq = {s_q:.3f}, sγ = {s_gamma:.3f}, sc = {s_c:.3f}")

    # Depth Factors (for D/B <= 1)
    d_q = 1 + 2 * tan_phi_d * (1 - math.sin(phi_d_prime_rad))**2 * (D / B)
    d_gamma = 1.0
    d_c = d_q - (1 - d_q) / (N_c * tan_phi_d)
    print(f"Depth Factors: dq = {d_q:.3f}, dγ = {d_gamma:.3f}, dc = {d_c:.3f}")

    # Inclination Factors
    m = (2 + B_L_ratio) / (1 + B_L_ratio)
    i_q = (1 - H_d / V_d)**m
    i_gamma = (1 - H_d / V_d)**(m + 1)
    i_c = (i_q * N_q - 1) / (N_q - 1)
    print(f"Inclination Factors: iq = {i_q:.3f}, iγ = {i_gamma:.3f}, ic = {i_c:.3f}")
    print("-" * 30)

    print("Step 5: Calculate Design Bearing Resistance (Rd)")
    # Overburden pressure
    q_prime = gamma_soil * D
    print(f"Effective Overburden Pressure (q') = {gamma_soil:.1f} kN/m³ * {D:.1f} m = {q_prime:.1f} kPa")

    # Individual terms of the bearing capacity equation
    term_c = c_d_prime * N_c * s_c * d_c * i_c
    term_q = q_prime * N_q * s_q * d_q * i_q
    term_gamma = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * d_gamma * i_gamma

    print("\nFinal Equation: R_d = A' * (c'd*Nc*sc*dc*ic + q'*Nq*sq*dq*iq + 0.5*γ*B'*Nγ*sγ*dγ*iγ)")
    print("Plugging in the values:")
    
    print(f"Term 1 (Cohesion): {c_d_prime:.2f} * {N_c:.2f} * {s_c:.3f} * {d_c:.3f} * {i_c:.3f} = {term_c:.2f} kPa")
    print(f"Term 2 (Surcharge): {q_prime:.1f} * {N_q:.2f} * {s_q:.3f} * {d_q:.3f} * {i_q:.3f} = {term_q:.2f} kPa")
    print(f"Term 3 (Self-weight): 0.5 * {gamma_soil:.1f} * {B_prime:.4f} * {N_gamma:.2f} * {s_gamma:.3f} * {d_gamma:.1f} * {i_gamma:.3f} = {term_gamma:.2f} kPa")

    # Total design resistance
    q_ult = term_c + term_q + term_gamma
    R_d = A_prime * q_ult
    
    print(f"\nTotal Bearing Pressure Resistance = {term_c:.2f} + {term_q:.2f} + {term_gamma:.2f} = {q_ult:.2f} kPa")
    print(f"Total Design Resistance (Rd) = {A_prime:.4f} m² * {q_ult:.2f} kPa = {R_d:.2f} kN")
    
    return R_d

# Run the calculation and print the final result
final_resistance = calculate_design_resistance()
print(f"\nFinal Answer: The design resistance is {final_resistance:.1f} kN.")
print(f"<<<{final_resistance:.1f}>>>")
