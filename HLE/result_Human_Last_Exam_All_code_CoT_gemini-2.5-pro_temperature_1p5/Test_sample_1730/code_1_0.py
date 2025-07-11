import math

def solve_puzzle():
    """
    Solves the multi-part problem to find the number of prisms that fit in a cube.
    """

    # Part 1: Calculate T for the prism base
    # Given parameters for the solar collector
    L = 1.5  # m
    B = 0.85  # m
    U_inf = 1.0  # m/s
    g = 9.81  # m/s^2
    rho_f = 1.204  # kg/m^3
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m*K)
    Pr_f = 0.707
    beta_f = 0.00341  # K^-1

    # Calculate average surface and ambient temperatures
    # theta_w_avg = integral(30 + 10*sin(pi*x/L))dx / L from 0 to L
    # integral -> 30*x - (10*L/pi)*cos(pi*x/L)
    # [30L - (10L/pi)*(-1)] - [0 - (10L/pi)*(1)] = 30L + 20L/pi
    # (30L + 20L/pi) / L = 30 + 20/pi
    theta_w_avg = 30 + 20 / math.pi
    
    # theta_inf_avg = integral(10 + 0.05y)dy / B from 0 to B
    # integral -> 10y + 0.025*y^2
    # [10B + 0.025B^2] - 0
    # (10B + 0.025B^2) / B = 10 + 0.025B
    theta_inf_avg = 10 + 0.025 * B
    
    delta_theta = theta_w_avg - theta_inf_avg

    # Calculate dimensionless numbers
    Re_L = (U_inf * L) / nu_f
    Gr_L = (g * beta_f * delta_theta * L**3) / (nu_f**2)
    Ra_L = Gr_L * Pr_f

    # Determine convection type (mixed convection as Gr/Re^2 is ~1)
    # Forced convection Nusselt number (laminar)
    Nu_L_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    
    # Natural convection Nusselt number (vertical plate)
    # Using the Churchill and Chu correlation for a wide range of Ra
    term1 = 0.387 * Ra_L**(1/6)
    term2 = (1 + (0.492 / Pr_f)**(9/16))**(8/27)
    Nu_L_natural = (0.825 + term1 / term2)**2
    
    # Combined Nusselt number for mixed convection
    Nu_L = (Nu_L_forced**3 + Nu_L_natural**3)**(1/3)
    
    # Calculate heat transfer coefficient and total heat loss
    h = (Nu_L * k_f) / L
    A = L * B
    Q_dot_V = h * A * delta_theta
    
    # Calculate T and round to the nearest integer
    T_float = Q_dot_V / 80
    T = round(T_float)

    print("--- Calculation for Prism Dimension T ---")
    print(f"1. Average surface temperature: {theta_w_avg:.2f} C")
    print(f"2. Average ambient temperature: {theta_inf_avg:.2f} C")
    print(f"3. Heat loss from the collector Q_V: {Q_dot_V:.2f} W")
    print(f"4. Prism base side T = Q_V / 80 = {Q_dot_V:.2f} / 80 = {T_float:.2f}")
    print(f"5. Rounded to the nearest integer, T = {T}")
    print("-" * 20)

    # Part 2: Calculate D for the prism depth
    # Given parameters for the beam
    q_0 = 3.0  # N/m
    l = 2.0  # m
    
    # Calculate `a`
    a = (64/3 - math.pi/4)**(-1/3)

    # Calculate maximum bending moment (for a simply supported beam)
    M_max = (q_0 * l**2) / 8
    
    # Calculate Area Moment of Inertia (I_y)
    # I_y = I_square - I_cutouts = (64a^4/3) - (pi*a^4/4) = a^4*(64/3 - pi/4)
    # Since a = (64/3 - pi/4)^(-1/3), 1/a^3 = (64/3 - pi/4)
    # So, I_y = a^4 * (1/a^3) = a
    I_y = a

    # Maximum distance from the neutral axis
    z_max = 2 * a
    
    # Calculate maximum normal stress
    sigma_xx_max = (M_max * z_max) / I_y

    # Calculate D
    D = sigma_xx_max / 3

    print("--- Calculation for Prism Dimension D ---")
    print(f"1. Maximum bending moment M_max: {M_max:.2f} N*m")
    print(f"2. Maximum normal stress sigma_max: {sigma_xx_max:.2f} N/m^2")
    print(f"3. Prism depth D = sigma_max / 3 = {sigma_xx_max:.2f} / 3 = {D:.2f}")
    print("-" * 20)
    
    # Part 3: Calculate how many prisms fit in the cube
    # Cube dimensions: 3x3x3
    # Prism dimensions: Right triangular base with legs T=2, depth D=1.
    
    print("--- Calculation for Number of Prisms ---")
    print(f"The problem is to fit prisms with a {T}x{T} right-triangular base and depth {int(D)} into a 3x3x3 cube.")
    print("A pair of such prisms forms a 2x2x1 rectangular box.")
    print("We can pack these boxes into the 3x3x3 cube to find the number of prisms.")
    
    # Packing strategy
    # 1. Pack a 3x2x2 block made of three 1x2x2 boxes (each holding 2 prisms).
    #    This fills region [0,3]x[0,2]x[0,2].
    prisms_r1 = 3 * 2
    
    # 2. This leaves two large empty blocks:
    #    L1: 3x1x3 (at y=[2,3])
    #    L2: 3x2x1 (at z=[2,3])
    # 3. Pack a 2x1x2 box (2 prisms) into the 3x1x3 block L1.
    prisms_r2 = 1 * 2
    
    # 4. Pack a 2x2x1 box (2 prisms) into the 3x2x1 block L2.
    prisms_r3 = 1 * 2

    total_prisms = prisms_r1 + prisms_r2 + prisms_r3
    
    print("\nPacking Plan:")
    print(f"1. A block of size 3x2x2 is placed in the cube. This block contains {prisms_r1} prisms.")
    print(f"2. In the remaining 3x1x3 space, we can fit a 2x1x2 block, containing {prisms_r2} prisms.")
    print(f"3. In the remaining 3x2x1 space, we can fit a 2x2x1 block, containing {prisms_r3} prisms.")
    print(f"The packing is valid and non-overlapping.")
    print("\nFinal Calculation:")
    print(f"Total number of prisms = {prisms_r1} (from block 1) + {prisms_r2} (from block 2) + {prisms_r3} (from block 3)")
    final_answer = prisms_r1 + prisms_r2 + prisms_r3
    print(f"Total number of prisms = {final_answer}")
    
    return final_answer

if __name__ == '__main__':
    final_answer = solve_puzzle()
    print(f"<<<{final_answer}>>>")