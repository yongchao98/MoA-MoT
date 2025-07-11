import math
from scipy.optimize import brentq

def calculate_and_print_solution():
    """
    Solves the footing design problem and prints the detailed solution.
    """
    # Step 1: Define given parameters
    Gk = 1000  # kN, characteristic permanent vertical load
    Qk_v = 1500  # kN, characteristic variable vertical load
    Qk_h = 300   # kN, characteristic variable horizontal load
    h_load_arm = 2 + 0.75 # m, height of horizontal load above footing base
    Df = 0.75    # m, depth of footing base
    gamma_concrete = 24 # kN/m^3
    gamma_soil = 20     # kN/m^3
    phi_k_deg = 35    # degrees, characteristic angle of shearing resistance
    c_prime_k = 0       # kPa, characteristic effective cohesion

    # ULS Partial Factors for DA1-C1
    gamma_G = 1.35 # permanent actions
    gamma_Q = 1.50 # variable actions
    gamma_phi = 1.0 # material property phi
    gamma_c = 1.0   # material property c
    gamma_gamma = 1.0 # material property unit weight

    # Design soil parameters
    phi_d_rad = math.atan(math.tan(math.radians(phi_k_deg)) / gamma_phi)
    phi_d_deg = math.degrees(phi_d_rad)
    c_prime_d = c_prime_k / gamma_c
    gamma_d = gamma_soil / gamma_gamma

    # Step 2: Define a function to find the equilibrium point (q_Rd - q_Ed = 0)
    def find_equilibrium_size(B):
        # Calculate Applied Loads and Stresses (q_Ed)
        W_base_k = B * B * Df * gamma_concrete
        Vd = gamma_G * (Gk + W_base_k) + gamma_Q * Qk_v
        Hd = gamma_Q * Qk_h
        Md = Hd * h_load_arm
        
        if Vd <= 0: return 1e6 # avoid division by zero
        e = Md / Vd
        
        if e >= B / 2: return 1e6 # Footing too small, no contact
            
        B_prime = B - 2 * e
        L_prime = B
        A_prime = B_prime * L_prime
        q_Ed = Vd / A_prime

        # Calculate Bearing Resistance (q_Rd)
        Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
        Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
        q_prime = gamma_d * Df
        sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
        sgamma = 1 - 0.3 * (B_prime / L_prime)
        dgamma = 1.0
        dq = 1 + 2 * math.tan(phi_d_rad) * (1 - math.sin(phi_d_rad))**2 * (Df / B)
        m = (2 + (B_prime/L_prime)) / (1 + (B_prime/L_prime))
        load_ratio = Hd / Vd
        if load_ratio >= 1: return 1e6

        iq = (1 - load_ratio)**m
        igamma = (1 - load_ratio)**(m + 1)
        
        term_q = q_prime * Nq * sq * dq * iq
        term_gamma = 0.5 * gamma_d * B_prime * Ngamma * sgamma * dgamma * igamma
        q_Rd = term_q + term_gamma
        
        return q_Rd - q_Ed

    # Step 3: Solve for B using a numerical solver
    try:
        B_optimal = brentq(find_equilibrium_size, 2.0, 3.0, xtol=1e-4)
    except ValueError:
        print("Solver failed to find a root. Please check the problem parameters.")
        return

    # Step 4: Calculate final values with the optimal B and print the results
    W_base_k = B_optimal**2 * Df * gamma_concrete
    Vd = gamma_G * (Gk + W_base_k) + gamma_Q * Qk_v
    Hd = gamma_Q * Qk_h
    Md = Hd * h_load_arm
    e = Md / Vd
    B_prime = B_optimal - 2 * e
    L_prime = B_optimal
    A_prime = B_prime * L_prime
    q_Ed = Vd / A_prime

    Nq = math.exp(math.pi * math.tan(phi_d_rad)) * (math.tan(math.radians(45) + phi_d_rad / 2))**2
    Ngamma = 2 * (Nq - 1) * math.tan(phi_d_rad)
    q_prime = gamma_d * Df
    sq = 1 + (B_prime / L_prime) * math.sin(phi_d_rad)
    sgamma = 1 - 0.3 * (B_prime / L_prime)
    dgamma = 1.0
    dq = 1 + 2 * math.tan(phi_d_rad) * (1 - math.sin(phi_d_rad))**2 * (Df / B_optimal)
    m = (2 + (B_prime/L_prime)) / (1 + (B_prime/L_prime))
    load_ratio = Hd / Vd
    iq = (1 - load_ratio)**m
    igamma = (1 - load_ratio)**(m + 1)
    
    term_q = q_prime * Nq * sq * dq * iq
    term_gamma = 0.5 * gamma_d * B_prime * Ngamma * sgamma * dgamma * igamma
    q_Rd = term_q + term_gamma
    
    print("This is an iterative design problem to find the footing size B where applied stress equals bearing resistance.\n")
    print(f"Optimal footing width B found to be: {B_optimal:.3f} m")
    print("-" * 60)
    print("CALCULATION OF APPLIED ULS STRESS (q_Ed)")
    print(f"Design Vertical Load (Vd) = 1.35 * ({Gk:.0f} + {W_base_k:.1f}) + 1.50 * {Qk_v:.0f} = {Vd:.1f} kN")
    print(f"Design Horizontal Load (Hd) = 1.50 * {Qk_h:.0f} = {Hd:.1f} kN")
    print(f"Design Moment (Md) = {Hd:.1f} * {h_load_arm:.2f} = {Md:.1f} kNm")
    print(f"Eccentricity (e) = {Md:.1f} / {Vd:.1f} = {e:.4f} m")
    print(f"Effective width (B') = {B_optimal:.3f} - 2 * {e:.4f} = {B_prime:.3f} m")
    print(f"Applied ULS Stress (q_Ed) = Vd / (B' * L) = {Vd:.1f} / ({B_prime:.3f} * {L_prime:.3f}) = {q_Ed:.1f} kN/m^2\n")

    print("-" * 60)
    print("CALCULATION OF REQUIRED ULS DESIGN BEARING RESISTANCE (q_Rd)")
    print("q_Rd = (q' * Nq * sq * dq * iq) + (0.5 * γ' * B' * Nγ * sγ * dγ * iγ)\n")
    print("Component values:")
    print(f"  q'     = {q_prime:.1f} kPa")
    print(f"  Nq     = {Nq:.2f}")
    print(f"  Nγ     = {Ngamma:.2f}")
    print(f"  sq     = {sq:.3f}")
    print(f"  sγ     = {sgamma:.3f}")
    print(f"  dq     = {dq:.3f}")
    print(f"  dγ     = {dgamma:.1f}")
    print(f"  iq     = {iq:.3f}")
    print(f"  iγ     = {igamma:.3f}\n")

    print("Final equation with numerical values:")
    print(f"q_Rd = ({q_prime:.1f} * {Nq:.2f} * {sq:.3f} * {dq:.3f} * {iq:.3f}) + "
          f"(0.5 * {gamma_d:.1f} * {B_prime:.3f} * {Ngamma:.2f} * {sgamma:.3f} * {dgamma:.1f} * {igamma:.3f})")
    print(f"q_Rd = {term_q:.1f} + {term_gamma:.1f} = {q_Rd:.1f} kN/m^2\n")

    print("The required ULS design bearing resistance is the value where the footing is stable (q_Ed ≈ q_Rd).")
    print(f"Required ULS Design Bearing Resistance = {q_Rd:.1f} kN/m^2")

calculate_and_print_solution()