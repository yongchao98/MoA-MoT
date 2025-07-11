import math
from scipy.optimize import fsolve

def calculate_bearing_resistance():
    """
    Calculates the required ULS design bearing resistance for a pad footing.
    
    This function solves for the footing dimension B where the design action (q_Ed)
    equals the design resistance (q_Rd) and then reports this value.
    The calculations are based on Eurocode 7, Annex D formulation.
    """
    # 1. Given parameters and partial factors (ULS, DA1-1)
    Gk = 1000.0  # kN (Permanent vertical load)
    Qk_v = 1500.0 # kN (Variable vertical load)
    Qk_h = 300.0  # kN (Variable horizontal load)
    h_lever_arm = 2.0 + 0.75 # m (Lever arm for horizontal load from base)
    
    Df = 0.75   # m (Depth of footing base)
    phi_k = 35.0  # degrees (Characteristic friction angle)
    c_k = 0.0     # kPa (Characteristic cohesion)
    gamma_soil = 20.0 # kN/m^3 (Unit weight of soil)
    gamma_conc = 24.0 # kN/m^3 (Unit weight of concrete)

    # Partial factors for actions (EC7 DA1, Combination 1)
    gamma_G = 1.35
    gamma_Q = 1.50
    
    # Partial factors for material properties (EC7 DA1, Set M1)
    gamma_phi = 1.0
    gamma_c = 1.0

    # 2. Design soil parameters
    phi_d_rad = math.atan(math.tan(math.radians(phi_k)) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)
    c_d = c_k / gamma_c
    
    # Horizontal Design Load (Hd)
    Hd = gamma_Q * Qk_h

    def objective_function(B_array):
        """
        Function to be solved: q_Ed - q_Rd = 0.
        B_array is a single-element array [B] from fsolve.
        """
        B = B_array[0]
        L = B  # Square footing

        # Calculate design loads including footing self-weight
        W_footing_k = (B * L * Df) * gamma_conc
        Vd = gamma_G * (Gk + W_footing_k) + gamma_Q * Qk_v
        Md = Hd * h_lever_arm
        
        # Effective width B'
        e = Md / Vd
        B_prime = B - 2 * e
        L_prime = L
        
        if B_prime <= 0:
            return 1e6 # Return a large number for invalid B

        A_prime = B_prime * L_prime
        
        # Bearing capacity factors (EC7 Annex D)
        Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45 + phi_d_deg / 2)))**2
        Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
        
        # Effective overburden pressure
        q_prime = gamma_soil * Df
        
        # Shape factors
        sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
        sgamma = 1 - 0.4 * (B_prime / L_prime)
        
        # Depth factors
        k = Df / B_prime
        if k > 1.0: # Hansen's modification for Df/B > 1
            k = math.atan(Df/B_prime)
        dq = 1 + 2 * math.tan(phi_d_rad) * (1 - math.sin(phi_d_rad))**2 * k
        dgamma = 1.0

        # Inclination factors
        m_exp = (2 + (B_prime/L_prime)) / (1 + (B_prime/L_prime))
        iq = (1 - Hd / Vd)**m_exp
        igamma = (1 - Hd / Vd)**(m_exp + 1)
        
        # Design Bearing Resistance (q_Rd)
        term_q = q_prime * Nq * sq * dq * iq
        term_gamma = 0.5 * gamma_soil * B_prime * Ngamma * sgamma * dgamma * igamma
        q_Rd = term_q + term_gamma

        # Design Bearing Pressure (q_Ed)
        q_Ed = Vd / A_prime
        
        return q_Ed - q_Rd

    # 3. Solve for B
    initial_guess = 2.5
    B_solution, = fsolve(objective_function, initial_guess)

    # 4. Calculate final values with the solved B
    B = B_solution
    L = B
    W_footing_k = (B * L * Df) * gamma_conc
    Vd = gamma_G * (Gk + W_footing_k) + gamma_Q * Qk_v
    Md = Hd * h_lever_arm
    e = Md / Vd
    B_prime = B - 2 * e
    L_prime = L
    A_prime = B_prime * L_prime
    
    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45 + phi_d_deg / 2)))**2
    Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
    q_prime = gamma_soil * Df
    sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
    sgamma = 1 - 0.4 * (B_prime / L_prime)
    k = Df / B_prime
    if k > 1.0: k = math.atan(Df/B_prime)
    dq = 1 + 2 * math.tan(phi_d_rad) * (1 - math.sin(phi_d_rad))**2 * k
    dgamma = 1.0
    m_exp = (2 + (B_prime/L_prime)) / (1 + (B_prime/L_prime))
    iq = (1 - Hd / Vd)**m_exp
    igamma = (1 - Hd / Vd)**(m_exp + 1)

    term_q = q_prime * Nq * sq * dq * iq
    term_gamma = 0.5 * gamma_soil * B_prime * Ngamma * sgamma * dgamma * igamma
    q_Rd = term_q + term_gamma

    # Print out the detailed calculation for the final resistance
    print("--- Geotechnical Design Calculation (ULS, Eurocode 7 DA1-1) ---\n")
    print("Objective: Find the required footing width 'B' where design action (q_Ed) equals design resistance (q_Rd).")
    print(f"Solved footing width B = {B:.2f} m\n")
    print("--- Calculation of ULS Design Bearing Resistance (q_Rd) for this footing ---\n")
    print(f"Design Vertical Load (Vd) = {Vd:.1f} kN")
    print(f"Design Horizontal Load (Hd) = {Hd:.1f} kN")
    print(f"Design Moment at base (Md) = {Md:.1f} kNm")
    print(f"Eccentricity (e) = {e:.3f} m")
    print(f"Effective Width (B') = {B_prime:.3f} m\n")

    print("The ULS design bearing resistance equation is:")
    print("q_Rd = (c'd * Nc * sc * dc * ic) + (q' * Nq * sq * dq * iq) + (0.5 * γ' * B' * Nγ * sγ * dγ * iγ)\n")
    print("Since cohesion c'd is 0, the first term is zero.")
    print("q_Rd = (q' * Nq * sq * dq * iq) + (0.5 * γ' * B' * Nγ * sγ * dγ * iγ)\n")
    print("Component values:")
    print(f"  q' (overburden) = γ_soil * Df = {gamma_soil} * {Df} = {q_prime:.1f} kPa")
    print(f"  γ' (soil unit weight) = {gamma_soil:.1f} kN/m^3")
    print(f"  B' (effective width) = {B_prime:.3f} m\n")
    
    print(f"  Nq (for φ'd={phi_d_deg:.1f}°) = {Nq:.2f}")
    print(f"  Nγ (for φ'd={phi_d_deg:.1f}°) = {Ngamma:.2f}\n")
    
    print(f"  sq (shape factor) = {sq:.3f}")
    print(f"  sγ (shape factor) = {sgamma:.3f}\n")

    print(f"  dq (depth factor) = {dq:.3f}")
    print(f"  dγ (depth factor) = {dgamma:.3f}\n")

    print(f"  iq (inclination factor) = {iq:.3f}")
    print(f"  iγ (inclination factor) = {igamma:.3f}\n")
    
    print("Final equation with values:")
    print(f"q_Rd = ({q_prime:.1f} * {Nq:.2f} * {sq:.3f} * {dq:.3f} * {iq:.3f}) + " \
          f"(0.5 * {gamma_soil:.1f} * {B_prime:.3f} * {Ngamma:.2f} * {sgamma:.3f} * {dgamma:.1f} * {igamma:.3f})")
    
    print(f"q_Rd = {term_q:.1f} + {term_gamma:.1f}")
    print(f"\nThe required ULS design bearing resistance is {q_Rd:.1f} kN/m^2")

if __name__ == '__main__':
    calculate_bearing_resistance()
