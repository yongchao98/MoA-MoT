import math

def calculate_bearing_resistance():
    """
    Calculates the ULS design bearing resistance for a pad footing based on Eurocode 7.
    An iterative process has determined that a footing of B=2.3m satisfies Vd = Rd.
    This script calculates the bearing resistance for B = 2.3m.
    """
    # 1. Input Parameters
    Gk = 1000  # kN (Characteristic permanent vertical load)
    Qk_v = 1500 # kN (Characteristic variable vertical load)
    Qk_h = 300  # kN (Characteristic variable horizontal load)
    h_load = 2.0  # m (Height of horizontal load application above ground)
    
    # Footing and Soil Properties
    B = 2.3     # m (Footing width, determined from iterative design, B=L)
    L = 2.3     # m (Footing length)
    D = 0.75    # m (Depth of footing base, also footing thickness)
    gamma_concrete = 24  # kN/m^3
    gamma_soil = 20      # kN/m^3
    phi_k_deg = 35       # degrees (Characteristic angle of shearing resistance)
    c_k = 0              # kPa (Characteristic effective cohesion)
    
    # Partial Factors for ULS Combination 1 (EN 1997-1, Set A1/M1)
    gamma_G = 1.35 # Partial factor for permanent actions
    gamma_Q = 1.50 # Partial factor for variable actions
    gamma_phi = 1.0 # Partial factor for tan(phi')
    gamma_c = 1.0   # Partial factor for c'
    
    print("Step 1: Calculate Design Actions (Vd, Hd, Md)")
    # Self-weight of footing base
    W_k = B * L * D * gamma_concrete
    # Total permanent load
    Gk_total = Gk + W_k
    # Design vertical load
    Vd = gamma_G * Gk_total + gamma_Q * Qk_v
    # Design horizontal load
    Hd = gamma_Q * Qk_h
    # Design moment at the footing base
    lever_arm = h_load + D
    Md = Hd * lever_arm
    
    print(f"Footing self-weight W_k = {B}m * {L}m * {D}m * {gamma_concrete}kN/m^3 = {W_k:.2f} kN")
    print(f"Design Vertical Load Vd = {gamma_G} * ({Gk} + {W_k:.2f}) + {gamma_Q} * {Qk_v} = {Vd:.2f} kN")
    print(f"Design Horizontal Load Hd = {gamma_Q} * {Qk_h} = {Hd:.2f} kN")
    print(f"Design Moment Md = {Hd:.2f} kN * ({h_load}m + {D}m) = {Md:.2f} kNm\n")

    print("Step 2: Calculate Effective Footing Dimensions")
    # Eccentricity
    e = Md / Vd
    # Effective width (moment is parallel to side B)
    B_prime = B - 2 * e
    L_prime = L
    A_prime = B_prime * L_prime
    print(f"Eccentricity e = {Md:.2f} kNm / {Vd:.2f} kN = {e:.4f} m")
    print(f"Effective width B' = {B}m - 2 * {e:.4f}m = {B_prime:.4f} m")
    print(f"Effective area A' = {B_prime:.4f}m * {L_prime}m = {A_prime:.4f} m^2\n")

    print("Step 3: Calculate Design Soil Parameters")
    # Design angle of shearing resistance
    phi_d_rad = math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)
    # Design cohesion
    c_d = c_k / gamma_c
    print(f"Design friction angle phi'_d = arctan(tan({phi_k_deg}°)/{gamma_phi}) = {phi_d_deg:.2f}°")
    print(f"Design cohesion c'_d = {c_k}kPa / {gamma_c} = {c_d:.2f} kPa\n")

    print("Step 4: Calculate Bearing Capacity, Shape, and Inclination Factors")
    # Bearing capacity factors (Vesic, as per EN 1997-1 Annex D)
    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45 + phi_d_deg / 2)))**2
    Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
    Nc = (Nq - 1) / math.tan(phi_d_rad) if math.tan(phi_d_rad) > 0 else 0
    print(f"Nq = {Nq:.2f}, Ngamma = {Ngamma:.2f}, Nc = {Nc:.2f}")

    # Shape factors
    sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
    sgamma = 1 - 0.3 * (B_prime / L_prime)
    sc = (sq * Nq - 1) / (Nq - 1) if (Nq-1) > 0.001 else 1.0 # For c'=0 this is not used
    print(f"sq = {sq:.3f}, sgamma = {sgamma:.3f}")
    
    # Inclination factors
    # Exponent m for load inclined parallel to B
    m = (2 + (B_prime / L_prime)) / (1 + (B_prime / L_prime))
    # Inclination factors
    Vd_plus_cot_term = Vd + A_prime * c_d / (math.tan(phi_d_rad) if phi_d_rad > 0 else 1e9)
    
    if Hd > Vd_plus_cot_term:
        iq, igamma = 0, 0 # if H > V, resistance is zero
    else:
        iq = (1 - Hd / Vd_plus_cot_term)**m
        igamma = (1 - Hd / Vd_plus_cot_term)**(m + 1)
        
    ic = iq - (1-iq)/(Nc * math.tan(phi_d_rad)) if Nc > 0 and math.tan(phi_d_rad) > 0 else 0 # for c'=0 not used

    print(f"Exponent m = {m:.3f}")
    print(f"iq = {iq:.3f}, igamma = {igamma:.3f}\n")

    print("Step 5: Calculate ULS Design Bearing Resistance (q_Rd)")
    # Effective overburden pressure at footing base
    q_prime = gamma_soil * D
    print(f"Effective overburden pressure q' = {gamma_soil} kN/m^3 * {D} m = {q_prime:.2f} kPa")
    
    # Cohesion term (will be zero as c'=0)
    term_c = c_d * Nc * sc * ic
    # Overburden term
    term_q = q_prime * Nq * sq * iq
    # Self-weight term
    term_gamma = 0.5 * gamma_soil * B_prime * Ngamma * sgamma * igamma
    
    # Total design bearing resistance per unit area (q_Rd)
    q_Rd = term_c + term_q + term_gamma
    
    print("\nFinal Bearing Resistance Calculation:")
    print(f"q_Rd = (c'*Nc*sc*ic) + (q'*Nq*sq*iq) + (0.5*g'*B'*Ng*sg*ig)")
    print(f"q_Rd = {term_c:.1f} + ({q_prime:.1f}*{Nq:.2f}*{sq:.3f}*{iq:.3f}) + (0.5*{gamma_soil:.1f}*{B_prime:.3f}*{Ngamma:.2f}*{sgamma:.3f}*{igamma:.3f})")
    print(f"q_Rd = {term_c:.1f} + {term_q:.1f} + {term_gamma:.1f}")
    print(f"The required ULS design bearing resistance is: {q_Rd:.1f} kN/m^2\n")

    # Final Check (Optional)
    Rd = q_Rd * A_prime
    print("--- Design Check ---")
    print(f"Total Resistance Rd = q_Rd * A' = {q_Rd:.1f} * {A_prime:.3f} = {Rd:.1f} kN")
    print(f"Design Load Vd = {Vd:.1f} kN")
    print(f"Utilization = Vd / Rd = {Vd/Rd:.3f}")
    print("The utilization is ~1.0, so the assumed B=2.3m and the resulting q_Rd are correct.")
    
    return q_Rd

# Execute the calculation and print the final answer in the required format
final_answer = calculate_bearing_resistance()
# print(f"<<<{final_answer:.1f}>>>")