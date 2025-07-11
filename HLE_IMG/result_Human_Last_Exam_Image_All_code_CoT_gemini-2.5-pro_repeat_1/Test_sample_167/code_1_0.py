import math

def solve_bearing_capacity():
    """
    Calculates the required ULS design bearing resistance for a pad footing
    according to Eurocode 7, Design Approach 1, Combination 1.
    """
    # 1. Given Data and Partial Factors
    Gk = 1000  # kN
    Qk_vert = 1500  # kN
    Qk_horiz = 300  # kN
    h_load = 2.0  # m (height of load application above ground)
    D = 0.75  # m (depth of footing)
    gamma_soil = 20  # kN/m^3
    phi_k_deg = 35  # degrees
    c_k = 0  # kPa

    # Partial factors for EC7 DA1-1
    gamma_G = 1.35
    gamma_Q = 1.5
    gamma_phi = 1.0
    gamma_c = 1.0

    print("Step 1: Calculate Design Actions (ULS Combination 1)")
    V_d = gamma_G * Gk + gamma_Q * Qk_vert
    H_d = gamma_Q * Qk_horiz
    M_d = H_d * (h_load + D)
    e = M_d / V_d
    print(f"Design Vertical Load (Vd) = {gamma_G} * {Gk} + {gamma_Q} * {Qk_vert} = {V_d:.2f} kN")
    print(f"Design Horizontal Load (Hd) = {gamma_Q} * {Qk_horiz} = {H_d:.2f} kN")
    print(f"Design Moment (Md) = {H_d:.2f} * ({h_load} + {D}) = {M_d:.2f} kNm")
    print(f"Eccentricity (e) = {M_d:.2f} / {V_d:.2f} = {e:.4f} m\n")

    print("Step 2: Calculate Design Soil Parameters")
    phi_d_deg = math.degrees(math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi))
    phi_d_rad = math.radians(phi_d_deg)
    c_d = c_k / gamma_c
    q_prime = gamma_soil * D
    print(f"Design Friction Angle (φ'd) = arctan(tan({phi_k_deg}°)/{gamma_phi}) = {phi_d_deg:.2f}°")
    print(f"Design Cohesion (c'd) = {c_k} / {gamma_c} = {c_d:.2f} kPa")
    print(f"Effective Overburden Pressure (q') = {gamma_soil} * {D} = {q_prime:.2f} kPa\n")
    
    print("Step 3: Calculate Bearing Capacity Factors")
    N_q = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.pi / 4 + phi_d_rad / 2) ** 2)
    N_gamma = 2 * (N_q - 1) * math.tan(phi_d_rad)
    print(f"For φ'd = {phi_d_deg:.2f}°, Nq = {N_q:.3f}, Nγ = {N_gamma:.3f}\n")

    def calculate_Rd(B):
        """Calculates the design resistance Rd for a given footing width B."""
        if e >= B / 6:
            # This check is for uplift, which should be avoided.
            # We return a very small number to guide the solver.
            return 1e-6 

        B_prime = B - 2 * e
        L_prime = B  # Since B=L
        A_prime = B_prime * L_prime
        
        ratio_B_L = B_prime / L_prime
        
        # Shape factors
        s_q = 1 + ratio_B_L * math.sin(phi_d_rad)
        s_gamma = 1 - 0.3 * ratio_B_L

        # Inclination factors
        m = (2 + ratio_B_L) / (1 + ratio_B_L)
        H_d_V_d_ratio = H_d / V_d
        i_q = (1 - H_d_V_d_ratio)**m
        i_gamma = (1 - H_d_V_d_ratio)**(m + 1)
        
        # Design bearing capacity (pressure)
        q_term = q_prime * N_q * s_q * i_q
        gamma_term = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * i_gamma
        q_Rd = q_term + gamma_term

        # Total design resistance (force)
        R_d = A_prime * q_Rd
        return R_d

    print("Step 4: Find required footing width 'B' by solving Vd = Rd")
    # Bisection method to find B
    low_B, high_B = 1.0, 5.0
    for _ in range(100): # Iterate for precision
        mid_B = (low_B + high_B) / 2
        if calculate_Rd(mid_B) < V_d:
            low_B = mid_B
        else:
            high_B = mid_B
    B_final = (low_B + high_B) / 2
    print(f"Iteration finished. Required footing width B ≈ {B_final:.3f} m\n")

    print("Step 5: Calculate final bearing resistance with the determined B")
    # Recalculate all values with the final B for printing
    B_prime = B_final - 2 * e
    L_prime = B_final
    A_prime = B_prime * L_prime
    ratio_B_L = B_prime / L_prime
    s_q = 1 + ratio_B_L * math.sin(phi_d_rad)
    s_gamma = 1 - 0.3 * ratio_B_L
    m = (2 + ratio_B_L) / (1 + ratio_B_L)
    H_d_V_d_ratio = H_d / V_d
    i_q = (1 - H_d_V_d_ratio)**m
    i_gamma = (1 - H_d_V_d_ratio)**(m + 1)
    
    q_term = q_prime * N_q * s_q * i_q
    gamma_term = 0.5 * gamma_soil * B_prime * N_gamma * s_gamma * i_gamma
    q_Rd = q_term + gamma_term
    R_d_final = A_prime * q_Rd

    print(f"Final parameters for B = {B_final:.3f} m:")
    print(f"  Effective width B' = {B_final:.3f} - 2 * {e:.4f} = {B_prime:.3f} m")
    print(f"  Effective area A' = {B_prime:.3f} * {L_prime:.3f} = {A_prime:.3f} m²")
    print(f"  Shape factor sq = 1 + ({B_prime:.3f}/{L_prime:.3f})*sin({phi_d_deg:.2f}) = {s_q:.3f}")
    print(f"  Shape factor sγ = 1 - 0.3*({B_prime:.3f}/{L_prime:.3f}) = {s_gamma:.3f}")
    print(f"  Inclination factor iq = (1 - {H_d:.2f}/{V_d:.2f})^{m:.3f} = {i_q:.3f}")
    print(f"  Inclination factor iγ = (1 - {H_d:.2f}/{V_d:.2f})^{m+1:.3f} = {i_gamma:.3f}\n")
    
    print("The ULS design bearing resistance (q_Rd) is given by the formula:")
    print("q_Rd = c'*Nc*sc*ic + q'*Nq*sq*iq + 0.5*γ'*B'*Nγ*sγ*iγ")
    print("Since c'=0, the first term is zero.")
    print(f"q_Rd = ({q_prime:.1f} * {N_q:.1f} * {s_q:.3f} * {i_q:.3f}) + (0.5 * {gamma_soil:.1f} * {B_prime:.3f} * {N_gamma:.1f} * {s_gamma:.3f} * {i_gamma:.3f})")
    print(f"q_Rd = {q_term:.1f} + {gamma_term:.1f}")
    print(f"q_Rd = {q_Rd:.1f} kN/m²\n")
    print("Check: Total resistance Rd = A' * q_Rd")
    print(f"Rd = {A_prime:.3f} m² * {q_Rd:.1f} kN/m² = {R_d_final:.1f} kN")
    print(f"This should be equal to the design load Vd = {V_d:.1f} kN. The small difference is due to rounding.")
    
    return q_Rd

if __name__ == '__main__':
    result = solve_bearing_capacity()
    print(f"\nFinal Answer: The required ULS design bearing resistance is {result:.1f} kN/m².")
    # The final answer block as requested
    print(f"\n<<<{result:.1f}>>>")
