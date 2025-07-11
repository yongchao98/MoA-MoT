import math

def solve_problem():
    """
    This function solves the entire multi-step problem.
    """
    # Part 1: Calculate T
    print("--- Part 1: Calculating T ---")
    
    # Given parameters for heat transfer
    L = 1.5  # m
    B = 0.85  # m
    U_inf = 1.0  # m/s
    rho_f = 1.204  # kg/m^3
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m.K)
    Pr_f = 0.707
    beta_f = 0.00341  # K^-1
    g = 9.81  # m/s^2

    # 1. Calculate average temperatures
    # theta_w_avg = (1/L) * integral(30 + 10*sin(pi*x/L)) dx from 0 to L
    # integral = [30x - (10L/pi)*cos(pi*x/L)] from 0 to L
    # = (30L + 10L/pi) - (-10L/pi) = 30L + 20L/pi
    # avg = 30 + 20/pi
    theta_w_avg = 30 + 20 / math.pi
    
    # theta_inf_avg = (1/B) * integral(10 + 0.05y) dy from 0 to B
    # integral = [10y + 0.025y^2] from 0 to B = 10B + 0.025B^2
    # avg = 10 + 0.025B
    theta_inf_avg = 10 + 0.025 * B
    
    delta_theta_avg = theta_w_avg - theta_inf_avg
    
    print(f"Average surface temperature (theta_w_avg): {theta_w_avg:.2f} C")
    print(f"Average ambient temperature (theta_inf_avg): {theta_inf_avg:.2f} C")
    print(f"Average temperature difference (delta_theta_avg): {delta_theta_avg:.2f} K")

    # 2. Determine convection regime
    # Forced convection along L, Natural convection on vertical height B
    Re_L = U_inf * L / nu_f
    Gr_B = (g * beta_f * delta_theta_avg * B**3) / (nu_f**2)
    
    print(f"Reynolds number (Re_L): {Re_L:.2e}")
    print(f"Grashof number (Gr_B): {Gr_B:.2e}")
    print(f"Gr/Re^2 ratio: {Gr_B / (Re_L**2):.2f}. This indicates mixed convection.")

    # 3. Calculate heat transfer coefficient
    # Forced convection (laminar flow correlation for flat plate)
    Nu_L_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    h_forced = (Nu_L_forced * k_f) / L
    
    # Natural convection (turbulent flow correlation for vertical plate, Ra > 10^9)
    Ra_B = Gr_B * Pr_f
    Nu_B_natural = 0.1 * (Ra_B**(1/3))
    h_natural = (Nu_B_natural * k_f) / B
    
    # Combine for mixed convection (perpendicular flows)
    h_mixed = (h_forced**3 + h_natural**3)**(1/3)
    
    print(f"Forced heat transfer coefficient (h_forced): {h_forced:.2f} W/(m^2.K)")
    print(f"Natural heat transfer coefficient (h_natural): {h_natural:.2f} W/(m^2.K)")
    print(f"Mixed heat transfer coefficient (h_mixed): {h_mixed:.2f} W/(m^2.K)")

    # 4. Calculate heat loss
    A = L * B
    Q_v = h_mixed * A * delta_theta_avg
    print(f"Total heat loss (Q_v): {Q_v:.2f} W")

    # 5. Calculate T
    T_float = Q_v / 80.0
    T = int(round(T_float))
    print(f"Calculated T = {Q_v:.2f} / 80 = {T_float:.2f}, which rounds to {T}")
    
    print("\n" + "="*40 + "\n")
    
    # Part 2: Calculate D
    print("--- Part 2: Calculating D ---")
    
    # Given parameters for beam problem
    q0 = 3.0  # N/m
    l_beam = 2.0  # m
    # a is given by a formula, but we find I_yy = a, so we don't need to calculate a's value.
    # a_val = (64/3 - math.pi/4)**(-1/3)
    
    # 1. Calculate max bending moment (simply supported beam)
    M_max = (q0 * l_beam**2) / 8
    print(f"Maximum bending moment (M_max): {M_max:.2f} N.m")
    
    # 2. Moment of inertia I_yy
    # As derived in the thought process, the complex integral for I_yy simplifies
    # I_yy = a^4 * (64/3 - pi/4).
    # And a = (64/3 - pi/4)^(-1/3), so a^3 = 1 / (64/3 - pi/4).
    # Therefore, I_yy = a * a^3 * (64/3 - pi/4) = a * 1 = a.
    # So, I_yy = a (in m^4, where a is in m).
    
    # 3. Calculate max normal stress
    # sigma_max = (M_max * z_max) / I_yy
    # z_max is the max distance from neutral axis, which is 2*a.
    # sigma_max = (M_max * 2*a) / a = 2 * M_max
    sigma_xx_max = 2 * M_max
    print(f"Maximum normal stress (sigma_xx_max): {sigma_xx_max:.2f} N/m^2")
    
    # 4. Calculate D
    D_float = sigma_xx_max / 3.0
    D = int(round(D_float))
    print(f"Calculated D = {sigma_xx_max:.2f} / 3 = {D_float:.2f}, which rounds to {D}")

    print("\n" + "="*40 + "\n")

    # Part 3: Calculate the number of prisms
    print("--- Part 3: Calculating the Number of Prisms ---")
    
    cube_dim = 3
    prism_leg = T
    prism_depth = D
    
    print(f"We need to fit prisms into a {cube_dim}x{cube_dim}x{cube_dim} cube.")
    print(f"Each prism has a right triangular base with legs of length {prism_leg}, and a depth of {prism_depth}.")
    
    # Since prism_depth is 1 and cube_dim is 3, we can stack 3 layers.
    num_layers = cube_dim // prism_depth
    print(f"The number of layers we can stack is {cube_dim} / {prism_depth} = {num_layers}.")
    
    # The problem reduces to packing 2x2 right triangles into a 3x3 square.
    # A valid, non-overlapping packing of 3 such triangles exists:
    # T1: vertices (0,0), (2,0), (0,2)
    # T2: vertices (3,3), (1,3), (3,1)
    # T3: vertices (1,1), (3,1), (1,3)
    prisms_per_layer = 3
    print(f"For the 2D packing problem, we can fit {prisms_per_layer} triangular bases into a {cube_dim}x{cube_dim} square.")
    
    total_prisms = prisms_per_layer * num_layers
    
    print("\n--- Final Calculation ---")
    print(f"Number of prisms = (Prisms per layer) * (Number of layers)")
    print(f"Number of prisms = {prisms_per_layer} * {num_layers} = {total_prisms}")
    
    print("\n" + "="*40 + "\n")
    print("Final Answer:")
    print(f"{total_prisms}")
    
    return total_prisms

# Execute the function to get the final answer
final_answer = solve_problem()
# The final answer is requested in a specific format.
# The code above already prints the final number.
# This is just to conform to the final output format requirement.
# print(f"<<<{final_answer}>>>")