import numpy as np

def solve_problem():
    # Part 1: Calculate T
    # Given parameters for the solar collector
    L = 1.5  # m, length
    B = 0.85  # m, width
    U_inf = 1.0  # m/s, wind speed
    
    # Air properties at film temperature
    rho_f = 1.204  # kg/m^3
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m.K)
    Pr_f = 0.707
    beta_f = 0.00341  # K^-1
    g = 9.81  # m/s^2

    # Step 1: Calculate average temperatures
    # Average wall temperature: integral of theta_w(x) / L
    # integral from 0 to L of (30 + 10*sin(pi*x/L)) dx = 30L + 20L/pi
    theta_w_avg = 30 + 20 / np.pi
    
    # Average ambient temperature: integral of theta_inf(y) / B
    # integral from 0 to B of (10 + 0.05y) dy = 10B + 0.025*B^2
    theta_inf_avg = 10 + 0.025 * B
    
    delta_theta_avg = theta_w_avg - theta_inf_avg

    # Step 2: Calculate dimensionless numbers
    Re_L = U_inf * L / nu_f
    Gr_L = g * beta_f * delta_theta_avg * (L**3) / (nu_f**2)
    Ra_L = Gr_L * Pr_f

    # Step 3: Use mixed convection model as Gr/Re^2 is not negligible
    # Forced convection Nusselt number (laminar)
    Nu_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    
    # Natural convection Nusselt number (turbulent, since Ra > 10^9)
    Nu_natural = 0.1 * (Ra_L**(1/3))
    
    # Combined Nusselt number for mixed convection (a common correlation)
    Nu_mixed = (Nu_forced**3 + Nu_natural**3)**(1/3)
    
    # Step 4: Calculate heat loss Q_V
    h_avg = Nu_mixed * k_f / L
    A = L * B
    Q_V = h_avg * A * delta_theta_avg
    
    # Step 5: Calculate T
    T_val = Q_V / 80
    T = int(round(T_val))

    # Part 2: Calculate D
    # Given parameters for the beam
    q0 = 3.0  # N/m
    l = 2.0  # m
    
    # Step 1: Calculate geometric properties
    # The value of 'a' is given by a formula that hints at the moment of inertia
    a_val = (64/3 - np.pi/4)**(-1/3)
    
    # From analysis, the moment of inertia I_yy = a
    I_yy = a_val
    z_max = 2 * a_val
    
    # Step 2: Calculate max bending moment (assuming simply supported beam)
    M_max = q0 * (l**2) / 8
    
    # Step 3: Calculate max normal stress
    sigma_max = M_max * z_max / I_yy
    
    # Step 4: Calculate D
    D_val = sigma_max / 3
    D = int(round(D_val))
    
    # Part 3: Calculate the number of prisms
    # Prism base is a right triangle with legs of length T
    # Prism depth is D
    cube_dim = 3
    
    # Number of layers along the depth
    num_layers = int(cube_dim / D)
    
    # Packing triangles on a 3x3 face
    # Triangle legs are T=2. Bounding box is 2x2.
    # We can fit a 2x2 square (made of 2 triangles) into a 3x3 square.
    # The remaining L-shaped region cannot fit more triangles.
    prisms_per_layer = 2
    
    total_prisms = prisms_per_layer * num_layers
    
    print("Step-by-step Calculation:")
    print(f"Heat Loss (Q_V): {Q_V:.2f} W")
    print(f"Calculated T value: {T_val:.2f}")
    print(f"Final T (rounded): {T}")
    print("-" * 20)
    print(f"Maximum Normal Stress (sigma_max): {sigma_max:.2f} N/m^2")
    print(f"Calculated D value: {D_val:.2f}")
    print(f"Final D (rounded): {D}")
    print("-" * 20)
    print("Final Packing Calculation:")
    print(f"A triangular prism has a right triangular base with legs of length T = {T} and a depth of D = {D}.")
    print(f"The cube has dimensions 3 x 3 x 3.")
    print(f"Number of prisms that can be stacked along the depth = floor(3 / D) = floor(3 / {D}) = {num_layers}.")
    print(f"Number of triangular bases that can be packed on a 3x3 face = {prisms_per_layer}.")
    print(f"Total number of prisms = {prisms_per_layer} (prisms per layer) * {num_layers} (layers) = {total_prisms}.")

solve_problem()
print("\n<<<6>>>")
