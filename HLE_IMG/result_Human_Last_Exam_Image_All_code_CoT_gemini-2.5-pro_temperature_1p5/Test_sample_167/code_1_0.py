import math

def calculate_bearing_resistance():
    """
    Calculates the required ULS design bearing resistance for a pad footing
    based on Eurocode 7 principles.
    """
    # Step 1: Define constants and parameters from the problem statement
    Gk = 1000.0  # kN, Characteristic permanent vertical load
    Qk_v = 1500.0 # kN, Characteristic variable vertical load
    Qk_h = 300.0  # kN, Characteristic variable horizontal load
    h_load = 2.0  # m, Lever arm for horizontal load
    D = 0.75      # m, Foundation embedment depth
    gamma_soil = 20.0  # kN/m^3, Unit weight of soil
    phi_k_deg = 35.0    # degrees, Characteristic angle of shearing resistance
    c_k = 0.0         # kPa, Characteristic effective cohesion

    # Partial factors for EC7 Design Approach 1, Combination 1
    gamma_G = 1.35
    gamma_Q = 1.50
    gamma_phi = 1.25
    gamma_c = 1.25

    # Step 2: Calculate design actions
    Vd = gamma_G * Gk + gamma_Q * Qk_v
    Hd = gamma_Q * Qk_h
    Md = Hd * h_load
    
    print("1. Design Actions (ULS Combination 1):")
    print(f"   - Design Vertical Load, Vd = {gamma_G} * {Gk} + {gamma_Q} * {Qk_v} = {Vd:.2f} kN")
    print(f"   - Design Horizontal Load, Hd = {gamma_Q} * {Qk_h} = {Hd:.2f} kN")
    print(f"   - Design Moment, Md = {Hd:.2f} * {h_load} = {Md:.2f} kNm")

    # Step 3: Calculate load eccentricity
    e = Md / Vd
    print(f"   - Load Eccentricity, e = {Md:.2f} / {Vd:.2f} = {e:.3f} m")
    
    # Step 4: Calculate design soil parameters
    phi_k_rad = math.radians(phi_k_deg)
    phi_d_rad = math.atan(math.tan(phi_k_rad) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)
    c_d = c_k / gamma_c
    
    print("\n2. Design Soil Parameters:")
    print(f"   - Design Friction Angle, φ'_d = arctan(tan({phi_k_deg}°)/{gamma_phi}) = {phi_d_deg:.2f}°")
    print(f"   - Design Cohesion, c'_d = {c_k} / {gamma_c} = {c_d:.2f} kPa")
    
    # Step 5: Calculate bearing capacity factors (Vesic's formulas)
    tan_phi_d = math.tan(phi_d_rad)
    Nq = math.exp(math.pi * tan_phi_d) * math.tan(math.radians(45) + phi_d_rad / 2)**2
    Ngamma = 2 * (Nq - 1) * tan_phi_d
    
    print("\n3. Bearing Capacity Factors (for φ'_d = {:.2f}°):".format(phi_d_deg))
    print(f"   - Nq = {Nq:.2f}")
    print(f"   - Nγ = {Ngamma:.2f}")
    
    # Effective overburden pressure at foundation base
    q_prime = gamma_soil * D
    
    # Step 6: Iterative search for footing width B
    def get_resistance_minus_load(B):
        B_prime = B - 2 * e
        if B_prime <= 0: return -Vd
        L_prime = B 
        A_prime = B_prime * L_prime
        
        # Shape factors (Vesic)
        sq = 1 + (B_prime / L_prime) * tan_phi_d
        sgamma = 1 - 0.4 * (B_prime / L_prime)
        
        # Depth factors (Hansen/Vesic)
        k = D / B
        dq = 1 + 2 * tan_phi_d * (1 - math.sin(phi_d_rad))**2 * k
        dgamma = 1.0
        
        # Inclination factors (Vesic)
        m = (2 + B_prime / L_prime) / (1 + B_prime / L_prime)
        iq = (1 - Hd / (Vd + A_prime * c_d / tan_phi_d if tan_phi_d > 0 else float('inf')))**m
        igamma = (1 - Hd / Vd)**(m + 1)
        
        # Ultimate bearing capacity (q_ult)
        term_q = q_prime * Nq * sq * dq * iq
        term_gamma = 0.5 * gamma_soil * B_prime * Ngamma * sgamma * dgamma * igamma
        q_ult = term_q + term_gamma
        
        # Total design resistance (Rd)
        Rd = q_ult * A_prime
        return Rd - Vd

    # Simple numerical search for B
    B_final = 0.1
    while get_resistance_minus_load(B_final) < 0:
        B_final += 0.001

    print(f"\n4. Iterative Sizing of Footing:")
    print(f"   - Required footing width, B, found to be ~{B_final:.2f} m to satisfy Vd <= Rd.")

    # Step 7: Calculate final bearing resistance with the determined B
    B_prime = B_final - 2 * e
    L_prime = B_final
    
    # Recalculate factors for the final B
    sq = 1 + (B_prime / L_prime) * tan_phi_d
    sgamma = 1 - 0.4 * (B_prime / L_prime)
    k = D / B_final
    dq = 1 + 2 * tan_phi_d * (1 - math.sin(phi_d_rad))**2 * k
    dgamma = 1.0
    m = (2 + B_prime / L_prime) / (1 + B_prime / L_prime)
    iq = (1 - Hd / (Vd + (B_prime * L_prime) * c_d / tan_phi_d if tan_phi_d > 0 else float('inf')))**m
    igamma = (1 - Hd / Vd)**(m + 1)
        
    term_q = q_prime * Nq * sq * dq * iq
    term_gamma = 0.5 * gamma_soil * B_prime * Ngamma * sgamma * dgamma * igamma
    required_resistance = term_q + term_gamma
    
    print("\n5. ULS Design Bearing Resistance Calculation (q_ult,d):")
    print("   The general bearing capacity equation (for c'=0) is:")
    print("   q_ult,d = (q' * Nq * sq * dq * iq) + (0.5 * γ' * B' * Nγ * sγ * dγ * iγ)\n")
    print("   With the following values:")
    print(f"   - Effective overburden pressure, q' = {gamma_soil} * {D} = {q_prime:.2f} kPa")
    print(f"   - Effective width, B' = {B_final:.2f} - 2 * {e:.3f} = {B_prime:.2f} m")
    print(f"   - Factors: Nq={Nq:.2f}, Nγ={Ngamma:.2f}, sq={sq:.3f}, sγ={sgamma:.3f}, dq={dq:.3f}, dγ={dgamma:.3f}, iq={iq:.3f}, iγ={igamma:.3f}")

    print("\n   Substituting these values into the two terms of the equation:")
    print(f"   Term 1 = ({q_prime:.2f} * {Nq:.2f} * {sq:.3f} * {dq:.3f} * {iq:.3f}) = {term_q:.2f} kN/m²")
    print(f"   Term 2 = (0.5 * {gamma_soil:.2f} * {B_prime:.2f} * {Ngamma:.2f} * {sgamma:.3f} * {dgamma:.3f} * {igamma:.3f}) = {term_gamma:.2f} kN/m²")
    
    print(f"\n   Required ULS design bearing resistance = {term_q:.2f} + {term_gamma:.2f} = {required_resistance:.2f} kN/m²")
    
    return required_resistance

if __name__ == '__main__':
    result = calculate_bearing_resistance()
    # The final answer part is for the platform, not for user copy-pasting.
    # The final print statement above gives the human-readable answer.
    # We round the final result to one decimal place for the answer block.
    # print(f"\nFinal Answer: {result:.1f}")
    # <<<526.9>>>

calculate_bearing_resistance()