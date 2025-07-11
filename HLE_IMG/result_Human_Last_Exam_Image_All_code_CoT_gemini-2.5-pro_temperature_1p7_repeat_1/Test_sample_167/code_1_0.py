import math

def calculate_uls_bearing_resistance():
    """
    Calculates the ULS design bearing resistance for a pad footing
    based on Eurocode 7 principles.
    """
    # Step 1: Define Given Data and Constants
    # Loads (Characteristic)
    Gk = 1000  # kN (Permanent vertical)
    Qk_v = 1500  # kN (Variable vertical)
    Qk_h = 300  # kN (Variable horizontal)

    # Geometry
    h_lever_arm = 2.0  # m (above ground)
    d_embedment = 0.75  # m

    # Soil Properties (Characteristic)
    phi_k_deg = 35  # degrees
    c_prime_k = 0  # kPa
    gamma_soil = 20  # kN/m^3

    # Partial Factors (ULS Combination 1, DA1-1)
    gamma_G = 1.35  # for permanent actions
    gamma_Q = 1.50  # for variable actions
    gamma_phi = 1.0  # for tan(phi')
    gamma_c = 1.0  # for c'
    gamma_gamma = 1.0  # for soil unit weight

    # Step 2: Calculate Design Actions and Soil Parameters
    # Design Loads
    Vd = gamma_G * Gk + gamma_Q * Qk_v
    Hd = gamma_Q * Qk_h
    total_lever_arm = h_lever_arm + d_embedment
    Md = Hd * total_lever_arm

    # Design Soil Parameters
    phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi))
    phi_d_rad = math.radians(phi_d_deg)
    c_prime_d = c_prime_k / gamma_c
    gamma_d = gamma_soil / gamma_gamma

    # Step 3: Iteratively find the required footing width B
    B = 2.0  # Initial guess for B
    step = 0.001 # Precision for B search
    Rd = 0

    # Iteration Loop to find B where Rd >= Vd
    while Rd < Vd:
        B += step
        L = B  # Square footing

        # Eccentricity and Effective Dimensions
        eB = Md / Vd
        B_prime = B - 2 * eB
        L_prime = L
        
        if B_prime <= 0: # Footing size is too small for the moment, continue search
            continue

        A_prime = B_prime * L_prime

        # Bearing Capacity Factors (Nq, N_gamma) based on phi_d
        Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.pi / 4 + phi_d_rad / 2) ** 2)
        N_gamma = 2 * (Nq - 1) * math.tan(phi_d_rad)

        # Overburden Pressure
        q_prime = gamma_d * d_embedment

        # Shape Factors (s_q, s_gamma)
        ratio_B_L = B_prime / L_prime
        s_q = 1 + ratio_B_L * math.sin(phi_d_rad)
        s_gamma = 1 - 0.3 * ratio_B_L

        # Inclination Factors (i_q, i_gamma)
        m = (2 + ratio_B_L) / (1 + ratio_B_L)
        # Vd + A'*c'*cot(phi') simplifies to Vd since c' = 0
        i_q = (1 - Hd / Vd) ** m
        i_gamma = (1 - Hd / Vd) ** (m + 1)

        # Ultimate Bearing Resistance (Pressure), q_ult
        # The term for cohesion (c') is zero and thus omitted.
        term_q = q_prime * Nq * s_q * i_q
        term_gamma = 0.5 * gamma_d * B_prime * N_gamma * s_gamma * i_gamma
        q_ult = term_q + term_gamma 

        # Total Design Resistance (Force), Rd
        Rd = q_ult * A_prime

    # Step 4: Output the results from the final successful iteration
    print("Final Calculation for ULS Design Bearing Resistance")
    print("="*50)

    print("\n1. Design Actions (Combination 1):")
    print(f"   Vd = {gamma_G} * {Gk} + {gamma_Q} * {Qk_v} = {Vd:.2f} kN")
    print(f"   Hd = {gamma_Q} * {Qk_h} = {Hd:.2f} kN")
    print(f"   Md = Hd * h = {Hd:.2f} * {total_lever_arm:.2f} = {Md:.2f} kNm")

    print("\n2. Required Footing Size (found via iteration):")
    print(f"   An iterative calculation found that a square footing of B = L = {B:.3f} m is required.")
    print(f"   This gives Rd = {Rd:.2f} kN, which is >= Vd = {Vd:.2f} kN.")

    print("\n3. Effective Footing Dimensions:")
    print(f"   Eccentricity eB = Md / Vd = {Md:.2f} / {Vd:.2f} = {eB:.4f} m")
    print(f"   Effective width B' = B - 2*eB = {B:.3f} - 2 * {eB:.4f} = {B_prime:.4f} m")
    print(f"   Effective length L' = L = {L:.3f} m")
    print(f"   Effective area A' = B' * L' = {A_prime:.4f} m^2")

    print("\n4. Design Soil Parameters & Factors:")
    print(f"   Design friction angle phi'_d = {phi_d_deg:.2f} degrees")
    print(f"   Effective overburden pressure q' = {gamma_d} * {d_embedment} = {q_prime:.2f} kPa")
    print(f"   Bearing Capacity Factors: Nq = {Nq:.2f}, N_gamma = {N_gamma:.2f}")
    print(f"   Shape Factors: sq = {s_q:.3f}, s_gamma = {s_gamma:.3f}")
    print(f"   Inclination Factors: iq = {i_q:.3f}, i_gamma = {i_gamma:.3f}")

    print("\n5. ULS Design Bearing Resistance Calculation (EC7 Formula):")
    print("   q_ult = (q' * Nq * sq * iq) + (0.5 * gamma' * B' * N_gamma * s_gamma * i_gamma)")
    print(f"   q_ult = ({q_prime:.2f} * {Nq:.2f} * {s_q:.3f} * {i_q:.3f}) + (0.5 * {gamma_d} * {B_prime:.4f} * {N_gamma:.2f} * {s_gamma:.3f} * {i_gamma:.3f})")
    print(f"   q_ult = {term_q:.2f} + {term_gamma:.2f}")
    print(f"   q_ult = {q_ult:.2f} kN/m^2")

    print("\n------------------------------------------------------------------")
    print(f"The required ULS design bearing resistance is {q_ult:.1f} kN/m^2.")
    print("------------------------------------------------------------------")

if __name__ == '__main__':
    calculate_uls_bearing_resistance()