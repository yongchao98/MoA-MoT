import math

def solve_bearing_resistance():
    """
    Calculates the required ULS design bearing resistance for a pad footing.
    This is an iterative process to find the footing size 'B' where the
    design bearing pressure (q_Ed) equals the design bearing resistance (q_Rd).
    """

    # 1. Define Input Parameters and Constants
    Gk = 1000.0  # Superimposed permanent vertical load (kN)
    Qk_v = 1500.0 # Variable vertical load (kN)
    Qk_h = 300.0  # Variable horizontal load (kN)
    
    h_col = 2.0   # Height of column/pier above footing (m)
    h_f = 0.75    # Footing thickness (m)
    Df = h_f      # Foundation depth (m)
    
    gamma_c = 24.0   # Unit weight of concrete (kN/m^3)
    gamma_soil = 20.0# Unit weight of soil (kN/m^3)
    
    phi_k_deg = 35.0 # Characteristic shear strength angle (degrees)
    c_prime_k = 0.0  # Characteristic effective cohesion (kPa)

    # Assumed column dimensions (not given in problem)
    b_col = 0.5   # Column width (m)
    l_col = 0.5   # Column length (m)

    # Eurocode 7 DA1-C1 partial factors
    gamma_G = 1.35  # Partial factor for permanent actions
    gamma_Q = 1.50  # Partial factor for variable actions
    gamma_phi = 1.0 # Partial factor for tan(phi')
    gamma_c_prime = 1.0 # Partial factor for c'
    gamma_gamma = 1.0 # Partial factor for soil unit weight
    
    # 2. Iterative Sizing to find optimal B
    
    # Design soil parameters
    phi_k_rad = math.radians(phi_k_deg)
    phi_d_rad = math.atan(math.tan(phi_k_rad) / gamma_phi)
    c_prime_d = c_prime_k / gamma_c_prime
    gamma_soil_d = gamma_soil / gamma_gamma
    
    # Bearing capacity factors (based on Vesic/EC7 formula)
    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
    # Ngamma for rough base
    Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
    Nc = (Nq - 1) / math.tan(phi_d_rad) if math.tan(phi_d_rad) > 0 else 0

    # Overburden pressure at formation level
    q_prime = gamma_soil_d * Df
    
    # Iteration using a bisection method to find B
    B_low = 0.1
    B_high = 5.0
    
    for _ in range(100): # Iterate to find the solution for B
        B = (B_low + B_high) / 2.0
        
        # 3. Calculate Design Actions
        W_pier = b_col * l_col * h_col * gamma_c
        W_footing = B * B * h_f * gamma_c
        Gk_total = Gk + W_pier + W_footing
        
        Vd = gamma_G * Gk_total + gamma_Q * Qk_v
        Hd = gamma_Q * Qk_h
        Md = Hd * (h_col + h_f)
        
        # 4. Determine Effective Area
        e = Md / Vd
        
        if e >= B / 6: # Check for tension, adjust search if needed
            B_low = B
            continue

        B_prime = B - 2 * e
        L_prime = B # Eccentricity in one direction
        A_prime = B_prime * L_prime
        
        # 5. Calculate Design Bearing Resistance (q_Rd)
        
        # Shape factors
        ratio_BL = B_prime / L_prime
        sq = 1 + ratio_BL * math.sin(phi_d_rad)
        sgamma = 1 - 0.3 * ratio_BL
        
        # Inclination factors
        # The term (Vd + A' * c'_d * cot(phi'_d)) simplifies to Vd since c'=0
        # Exponent 'm' for load inclination
        m = (2 + ratio_BL) / (1 + ratio_BL)
        
        iq = (1 - Hd / Vd)**m
        igamma = (1 - Hd / Vd)**(m + 1)
        
        # Calculate resistance q_Rd (c' term is zero)
        q_term = q_prime * Nq * sq * iq
        gamma_term = 0.5 * gamma_soil_d * B_prime * Ngamma * sgamma * igamma
        q_Rd = q_term + gamma_term

        # Calculate design bearing pressure (action)
        q_Ed = Vd / A_prime
        
        # Update search range
        if q_Ed > q_Rd:
            B_low = B
        else:
            B_high = B

    # 6. Final Calculation with the determined B
    B_final = B_high
    
    # Recalculate all final values with B_final
    W_pier = b_col * l_col * h_col * gamma_c
    W_footing = B_final * B_final * h_f * gamma_c
    Gk_total = Gk + W_pier + W_footing
    Vd = gamma_G * Gk_total + gamma_Q * Qk_v
    Hd = gamma_Q * Qk_h
    Md = Hd * (h_col + h_f)
    e = Md / Vd
    B_prime = B_final - 2 * e
    L_prime = B_final
    A_prime = B_prime * L_prime
    ratio_BL = B_prime / L_prime
    sq = 1 + ratio_BL * math.sin(phi_d_rad)
    sgamma = 1 - 0.3 * ratio_BL
    m = (2 + ratio_BL) / (1 + ratio_BL)
    iq = (1 - Hd / Vd)**m
    igamma = (1 - Hd / Vd)**(m + 1)
    q_term = q_prime * Nq * sq * iq
    gamma_term = 0.5 * gamma_soil_d * B_prime * Ngamma * sgamma * igamma
    required_resistance = q_term + gamma_term

    # 7. Print the final results and equation
    print("The required ULS design bearing resistance is found by sizing the footing such that the applied pressure equals the soil's resistance.")
    print(f"The optimal footing size was found to be B = {B_final:.2f} m.\n")
    print("The ULS design bearing resistance, q_R;d, is calculated using the formula:")
    print("q_R;d = (c' * Nc * sc * ic) + (q' * Nq * sq * iq) + (0.5 * γ' * B' * Nγ * sγ * iγ)\n")
    print(f"Since c' = 0, the first term is zero.")
    print("The final calculation is based on the following values:")
    print(f"  q' (overburden pressure) = {q_prime:.2f} kPa")
    print(f"  γ' (soil unit weight) = {gamma_soil_d:.2f} kN/m^3")
    print(f"  B' (effective width) = {B_prime:.2f} m")
    print(f"  Nq = {Nq:.2f}, Nγ = {Ngamma:.2f}")
    print(f"  sq = {sq:.3f}, sγ = {sgamma:.3f}")
    print(f"  iq = {iq:.3f}, iγ = {igamma:.3f}\n")
    
    print("Substituting these values into the second and third terms of the equation:")
    print(f"q' term = {q_prime:.2f} * {Nq:.2f} * {sq:.3f} * {iq:.3f} = {q_term:.2f} kPa")
    print(f"γ' term = 0.5 * {gamma_soil_d:.2f} * {B_prime:.2f} * {Ngamma:.2f} * {sgamma:.3f} * {igamma:.3f} = {gamma_term:.2f} kPa")
    
    print("\nRequired ULS Design Bearing Resistance (q_R;d) = q' term + γ' term")
    print(f"q_R;d = {q_term:.2f} + {gamma_term:.2f} = {required_resistance:.2f} kN/m^2")
    
    # Return the final value for the answer tag
    return required_resistance

if __name__ == '__main__':
    result = solve_bearing_resistance()
    # The final answer is formatted as requested at the end.
    # print(f"\n<<< {result:.1f} >>>")