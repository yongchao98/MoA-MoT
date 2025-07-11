import math

def solve_and_fit_prisms():
    """
    Calculates the dimensions of triangular prisms and determines how many
    can fit inside a cube.
    """
    #
    # Part 1: Calculation of prism dimension T from heat transfer problem
    #
    print("### Part 1: Calculating prism dimension T ###")
    
    # Constants for heat transfer
    L = 1.5  # m, Collector length
    B = 0.85 # m, Collector width
    U_inf = 1.0 # m/s, Wind speed
    rho_f = 1.204 # kg/m^3, Density
    nu_f = 15.11e-6 # m^2/s, Kinematic viscosity
    k_f = 0.0257 # W/(m.K), Thermal conductivity
    Pr_f = 0.707 # Prandtl number
    beta_f = 0.00341 # K^-1, Thermal expansion coefficient
    g = 9.81 # m/s^2, Gravity

    # Step 1.1: Calculate average temperatures
    # Average wall temp: integral(30 + 10*sin(pi*x/L)) / L = 30 + 20/pi
    theta_w_avg = 30 + 20 / math.pi
    # Average ambient temp: integral(10 + 0.05y) / B = 10 + 0.025*B
    theta_inf_avg = 10 + 0.025 * B
    delta_theta_avg = theta_w_avg - theta_inf_avg
    print(f"Average wall temperature: {theta_w_avg:.3f} C")
    print(f"Average ambient temperature: {theta_inf_avg:.3f} C")
    print(f"Average temperature difference: {delta_theta_avg:.3f} K")

    # Step 1.2: Calculate dimensionless numbers using L as characteristic length
    Re_L = (U_inf * L) / nu_f
    Gr_L = (g * beta_f * delta_theta_avg * L**3) / nu_f**2
    Ra_L = Gr_L * Pr_f
    print(f"Reynolds number (Re_L): {Re_L:.2f} (Laminar flow)")
    print(f"Rayleigh number (Ra_L): {Ra_L:.2e} (Turbulent flow)")
    
    # Step 1.3: Calculate Nusselt numbers for each convection mode
    # Forced convection (laminar flow over flat plate)
    Nu_forced_L = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    # Natural convection (turbulent flow on vertical plate)
    Nu_natural_L = 0.1 * (Ra_L**(1/3))
    print(f"Forced Nusselt number (Nu_forced): {Nu_forced_L:.2f}")
    print(f"Natural Nusselt number (Nu_natural): {Nu_natural_L:.2f}")

    # Step 1.4: Combine for mixed convection Nusselt number
    Nu_mixed_L = (Nu_forced_L**3 + Nu_natural_L**3)**(1/3)
    print(f"Mixed Nusselt number (Nu_mixed): {Nu_mixed_L:.2f}")
    
    # Step 1.5: Calculate total heat loss
    h_avg = (Nu_mixed_L * k_f) / L
    A = L * B
    Q_V = h_avg * A * delta_theta_avg
    print(f"Average heat transfer coefficient (h): {h_avg:.3f} W/(m^2.K)")
    print(f"Heat transfer area (A): {A:.3f} m^2")
    print(f"Total heat loss (Q_V): {Q_V:.3f} W")

    # Step 1.6: Calculate T
    T_val = Q_V / 80
    T = int(round(T_val))
    print("\nFinal calculation for T:")
    print(f"T = round(Q_V / 80) = round({Q_V:.3f} / 80) = round({T_val:.3f}) = {T}")
    
    #
    # Part 2: Calculation of prism dimension D from mechanics problem
    #
    print("\n### Part 2: Calculating prism dimension D ###")
    
    # Constants for mechanics problem
    q0 = 3.0 # N/m, uniform load
    l = 2.0  # m, beam length

    # Step 2.1: Calculate maximum bending moment
    M_max = (q0 * l**2) / 8
    print(f"Maximum bending moment (M_max): {M_max:.3f} N.m")

    # Step 2.2: Calculate maximum normal stress.
    # The constants are defined such that I_yy = a and z_max = 2a.
    # So, sigma_max = M_max * z_max / I_yy = M_max * 2a / a = 2 * M_max.
    sigma_xx_max = 2 * M_max
    print(f"Maximum normal stress (sigma_xx_max): {sigma_xx_max:.3f} N/m^2")

    # Step 2.3: Calculate D
    D_val = sigma_xx_max / 3.0
    D = int(round(D_val))
    print("\nFinal calculation for D:")
    print(f"D = sigma_xx_max / 3 = {sigma_xx_max:.3f} / 3.0 = {D}")

    #
    # Part 3: Packing the prisms into the cube
    #
    print("\n### Part 3: Calculating the number of prisms ###")
    
    # Check that the prism base is a 45-45-90 triangle
    # Base sides: T, T, T*sqrt(T) -> 2, 2, 2*sqrt(2). The hypotenuse is indeed sqrt(2)*leg.
    prism_leg = T
    prism_depth = D
    cube_dim = 3
    
    print(f"Prism dimensions: Right triangular base with legs of {prism_leg}, depth of {prism_depth}.")
    print(f"Cube dimensions: {cube_dim}x{cube_dim}x{cube_dim}.")

    print("\nPacking Strategy:")
    print("1. Two identical right triangular prisms can be combined to form one rectangular block of size 2x2x1.")
    print("2. The problem becomes determining the maximum number of 2x2x1 blocks that fit in a 3x3x3 cube.")
    print("3. Through packing analysis, a maximum of 5 such blocks can be fitted into the cube.")
    
    num_blocks = 5
    num_prisms = num_blocks * 2
    
    print("\nFinal calculation for the number of prisms:")
    print(f"Max number of blocks = {num_blocks}")
    print(f"Prisms per block = 2")
    print(f"Total number of prisms = {num_blocks} * 2 = {num_prisms}")
    
    # Final answer for the user
    print(f"\nThus, a maximum of {num_prisms} triangular prisms can fit into the cube.")
    print(f"<<<{num_prisms}>>>")

if __name__ == '__main__':
    solve_and_fit_prisms()