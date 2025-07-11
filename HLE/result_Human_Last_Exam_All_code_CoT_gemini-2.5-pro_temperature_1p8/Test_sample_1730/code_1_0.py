import math

def solve_puzzle():
    """
    This function solves the entire multi-step problem.
    1. Calculates the prism dimension T based on a heat transfer problem.
    2. Calculates the prism dimension D based on a beam mechanics problem.
    3. Calculates how many prisms fit into a 3x3x3 cube.
    """

    # Part 1: Calculate T from heat transfer
    
    # Given parameters
    L = 1.5         # m, plate length
    B = 0.85        # m, plate width (vertical height)
    U_inf = 1.0     # m/s, wind speed
    g = 9.81        # m/s^2

    # Fluid properties at film temperature
    k_f = 0.0257    # W/(m*K), thermal conductivity
    nu_f = 15.11e-6 # m^2/s, kinematic viscosity
    Pr_f = 0.707    # Prandtl number
    beta_f = 0.00341# K^-1, thermal expansion coefficient

    # Calculate average wall temperature by integrating theta_w(x)
    # integral[30 + 10*sin(pi*x/L)]dx from 0 to L is 30*L + 20*L/pi
    theta_w_avg = 30 + 20 / math.pi  # °C

    # Calculate average ambient temperature by integrating theta_inf(y)
    # integral[10 + 0.05*y]dy from 0 to B is 10*B + 0.025*B^2
    theta_inf_avg = 10 + 0.025 * B  # °C

    # Average temperature difference for convection
    delta_T_avg = theta_w_avg - theta_inf_avg # K

    # Forced convection calculations
    Re_L = U_inf * L / nu_f
    Nu_L_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3.0))
    h_forced = Nu_L_forced * k_f / L

    # Natural convection calculations
    Gr_B = (g * beta_f * delta_T_avg * (B**3)) / (nu_f**2)
    Ra_B = Gr_B * Pr_f
    # Correlation for turbulent natural convection on a vertical plate
    Nu_B_natural = 0.1 * (Ra_B**(1/3.0))
    h_natural = Nu_B_natural * k_f / B

    # Combined mixed convection coefficient using Churchill's correlation
    h_mixed = (h_forced**3 + h_natural**3)**(1/3.0)

    # Total heat loss
    Area = L * B
    Q_V = h_mixed * Area * delta_T_avg

    # Calculate T, rounded to the nearest integer
    T = int(round(Q_V / 80.0))

    # Part 2: Calculate D from beam mechanics
    
    # Given parameters
    q0 = 3.0  # N/m, distributed load
    l = 2.0   # m, beam length

    # Assuming a simply supported beam, the max moment is at the center
    M_max = q0 * (l**2) / 8.0

    # The problem setup for the cross-section (with a= (64/3-pi/4)^(-1/3)) 
    # leads to the numerical value of I_y being equal to a, while Z_max=2a.
    # sigma_xx_max = M_max * Z_max / I_y = M_max * 2a / a = 2 * M_max.
    sigma_xx_max = 2 * M_max

    # Calculate D
    D = int(round(sigma_xx_max / 3.0))

    # Part 3: Calculate how many prisms fit in the cube

    cube_dim = 3
    
    # A triangular prism with base T x T and depth D fits in a bounding box T x T x D.
    # We find how many such boxes fit in an axis-aligned configuration.
    num_boxes_x = math.floor(cube_dim / T)
    num_boxes_y = math.floor(cube_dim / T)
    num_boxes_z = math.floor(cube_dim / D)
    total_boxes = num_boxes_x * num_boxes_y * num_boxes_z
    
    # Two prisms fit perfectly into one bounding box.
    prisms_per_box = 2
    num_prisms = total_boxes * prisms_per_box

    # Print the results step-by-step
    print(f"Step 1: The calculated value for T is {T}.")
    print(f"Step 2: The calculated value for D is {D}.")
    print(f"Step 3: Calculate the number of triangular prisms that can fit in a {cube_dim}x{cube_dim}x{cube_dim} cube.")
    print(f"The prism's right triangular base has side lengths {T}, {T}, and {T}*sqrt({T}) ({T*math.sqrt(T):.2f}), with a depth of {D}.")
    print("To find the number of prisms, we calculate how many axis-aligned bounding boxes of size 2x2x1 fit into the 3x3x3 cube, and multiply by 2 (since 2 prisms fit per box).")
    print(f"\nFinal Calculation:")
    print(f"Number of Prisms = floor({cube_dim}/{T}) * floor({cube_dim}/{T}) * floor({cube_dim}/{D}) * {prisms_per_box}")
    print(f"= {num_boxes_x} * {num_boxes_y} * {num_boxes_z} * {prisms_per_box}")
    print(f"= {num_prisms}")

solve_puzzle()