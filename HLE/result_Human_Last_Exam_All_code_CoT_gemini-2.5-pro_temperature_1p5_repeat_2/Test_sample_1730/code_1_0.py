import math

def solve_problem():
    # Part 1: Calculate T
    
    # Given parameters for heat transfer
    L = 1.5  # m
    B = 0.85  # m
    U_inf = 1.0  # m/s
    g = 9.81  # m/s^2
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m.K)
    Pr_f = 0.707
    beta_f = 0.00341  # K^-1
    
    # Step 1.1 & 1.2: Calculate average temperatures
    # theta_w_avg = (1/L) * integral(30 + 10*sin(pi*x/L)) from 0 to L
    # = 30 + (10/L) * [-L/pi * cos(pi*x/L)] from 0 to L
    # = 30 + (10/pi) * (-cos(pi) - (-cos(0))) = 30 + 20/pi
    theta_w_avg = 30 + 20 / math.pi
    
    # theta_inf_avg = (1/B) * integral(10 + 0.05y) from 0 to B
    # = (1/B) * [10y + 0.025y^2] from 0 to B = 10 + 0.025*B
    theta_inf_avg = 10 + 0.025 * B
    
    delta_T_avg = theta_w_avg - theta_inf_avg
    
    # Step 1.3: Characterize flow
    Re_L = U_inf * L / nu_f
    Gr_L = g * beta_f * delta_T_avg * (L**3) / (nu_f**2)
    
    # Step 1.4: Calculate Nusselt number for mixed convection
    # Forced convection for laminar flow over a flat plate
    Nu_L_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    
    # Natural convection for a vertical plate
    Ra_L = Gr_L * Pr_f
    Nu_L_natural_term = 0.387 * (Ra_L**(1/6)) / ((1 + (0.492 / Pr_f)**(9/16))**(8/27))
    Nu_L_natural = (0.825 + Nu_L_natural_term)**2
    
    # Combined mixed convection Nusselt number
    Nu_L = (Nu_L_forced**3 + Nu_L_natural**3)**(1/3)
    
    # Step 1.5 & 1.6: Calculate heat loss
    h_avg = Nu_L * k_f / L
    Area = L * B
    Q_V = h_avg * Area * delta_T_avg
    
    # Step 1.7: Calculate T
    T_unrounded = Q_V / 80.0
    T = int(round(T_unrounded))
    
    # Part 2: Calculate D

    # Given parameters for beam problem
    q0 = 3.0  # N/m
    l = 2.0  # m
    
    # Step 2.1: Calculate maximum bending moment for a simply supported beam
    M_max = q0 * (l**2) / 8.0
    
    # Step 2.2 & 2.3: The problem is set up such that sigma_xx_max simplifies.
    # The value 'a' is defined as a = (64/3 - pi/4)^(-1/3).
    # The moment of inertia I_y for the described cross-section is I_y = (64/3 - pi/4)*a^4.
    # Substituting the definition of a, we get a^3 = 1/(64/3 - pi/4), or (64/3-pi/4)*a^3 = 1.
    # So, I_y = ((64/3 - pi/4)*a^3)*a = 1*a = a.
    # The maximum stress sigma = M_max * z_max / I_y, with z_max = 2a.
    # sigma_xx_max = M_max * (2a) / a = 2 * M_max.
    # Note: A full derivation shows units work out to N/m^2.
    # sigma_xx_max = (1.5 N.m) * (2*a m) / (a m^4) = 3 N/m^2.
    sigma_xx_max = 2 * M_max
    
    # Step 2.4: Calculate D
    D = sigma_xx_max / 3.0
    
    # Part 3: Calculate the number of prisms
    
    # Cube side length
    cube_side = 3.0
    
    # Prism base is a T x T right triangle, depth is D.
    # Bounding box of the prism is T x T x D.
    
    num_boxes_T1 = math.floor(cube_side / T)
    num_boxes_T2 = math.floor(cube_side / T)
    num_boxes_D = math.floor(cube_side / D)
    
    num_bounding_boxes = num_boxes_T1 * num_boxes_T2 * num_boxes_D
    
    # Two prisms fit into one bounding box
    num_prisms = 2 * num_bounding_boxes
    
    # Print the results
    print(f"Calculated T = {T_unrounded:.2f}, which is rounded to {T}")
    print(f"Calculated D = {D:.2f}")
    print("\nThe triangular prism has a right triangular base with perpendicular sides of length T={T} and a depth of D={D}.".format(T=T, D=D))
    print("The cube has dimensions 3 x 3 x 3.")
    
    print("\nTo find how many prisms fit, we calculate the number of bounding boxes that can be packed.")
    print("The number of prisms is calculated based on the formula: 2 * floor(cube_side / T) * floor(cube_side / T) * floor(cube_side / D)")
    
    final_eq_str = f"2 * floor({cube_side} / {T}) * floor({cube_side} / {T}) * floor({cube_side} / {D})"
    
    print(f"Number of prisms = {final_eq_str} = {num_prisms}")

solve_problem()
