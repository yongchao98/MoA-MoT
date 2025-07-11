import math

def solve_puzzle():
    """
    This function solves the entire multi-step problem by:
    1. Calculating the heat transfer to find the value of T.
    2. Calculating the beam stress to find the value of D.
    3. Solving the packing problem to find the number of prisms.
    """

    # Part 1: Calculate T from the heat transfer problem
    # Constants
    L = 1.5         # m, Length of the collector
    B = 0.85        # m, Width (height) of the collector
    U_inf = 1.0     # m/s, Wind speed
    g = 9.81        # m/s^2, Gravity
    nu_f = 15.11e-6 # m^2/s, Kinematic viscosity
    k_f = 0.0257    # W/(m.K), Thermal conductivity
    Pr_f = 0.707    # Prandtl number
    beta_f = 0.00341# K^-1, Thermal expansion coefficient

    # Calculate average surface and ambient temperatures
    # theta_w_avg = (1/L) * integral(30 + 10*sin(pi*x/L)) dx from 0 to L
    theta_w_avg = 30.0 + 10.0 * (2.0 / math.pi)
    # theta_inf_avg = (1/B) * integral(10 + 0.05*y) dy from 0 to B
    theta_inf_avg = 10.0 + 0.05 * B / 2.0
    
    delta_theta = theta_w_avg - theta_inf_avg

    # Calculate dimensionless numbers for mixed convection
    # Reynolds number (forced convection)
    Re_L = U_inf * L / nu_f
    # Forced convection Nusselt number (laminar flow over flat plate)
    Nu_F_L = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))

    # Rayleigh number (natural convection)
    Ra_L = (g * beta_f * delta_theta * (L**3) * Pr_f) / (nu_f**2)
    # Natural convection Nusselt number (Churchill and Chu correlation for vertical plate)
    term1_Nu_N = 0.387 * (Ra_L**(1/6))
    term2_Nu_N = (1 + (0.492 / Pr_f)**(9/16))**(8/27)
    Nu_N_L = (0.825 + term1_Nu_N / term2_Nu_N)**2

    # Combined Nusselt number for mixed convection (n=3)
    Nu_L = (Nu_F_L**3 + Nu_N_L**3)**(1/3)

    # Calculate heat transfer coefficient and total heat loss
    h = Nu_L * k_f / L
    A = L * B
    Q_V = h * A * delta_theta

    # Calculate T and round to the nearest integer
    T_float = Q_V / 80.0
    T = int(round(T_float))

    # Part 2: Calculate D from the beam mechanics problem
    # Constants
    q0 = 3.0     # N/m, Uniformly distributed load
    l_beam = 2.0 # m, Length of the beam

    # Based on the problem's formulation, the values are constructed to simplify perfectly.
    # M_max = q0 * l_beam^2 / 8 = 3 * 2^2 / 8 = 1.5 N.m
    # The term for 'a' implies I_y = a and z_max = 2a.
    # So, sigma_xx_max = (M_max * z_max) / I_y = (1.5 * 2*a) / a = 3.0 N/m^2.
    sigma_xx_max = 3.0
    
    # Calculate D
    D_float = sigma_xx_max / 3.0
    D = int(round(D_float))

    # Part 3: Solve the packing problem
    # Container is a cube with side length C
    C = 3
    
    # We form a rectangular block of size T x T x D by pairing two prisms.
    # For T=2, the base is a right triangle with legs T=2, so this is valid.
    
    # Find the maximum number of blocks that can fit by checking orientations.
    # Orientation 1: (T, T, D)
    n_blocks_1 = math.floor(C/T) * math.floor(C/T) * math.floor(C/D)
    # Orientation 2: (T, D, T)
    n_blocks_2 = math.floor(C/T) * math.floor(C/D) * math.floor(C/T)
    # Orientation 3: (D, T, T)
    n_blocks_3 = math.floor(C/D) * math.floor(C/T) * math.floor(C/T)
    
    num_blocks = max(n_blocks_1, n_blocks_2, n_blocks_3)

    # Total number of prisms is 2 * num_blocks
    num_prisms = 2 * num_blocks

    # Output the results and the final calculation as requested
    print(f"The calculated integer value for T is {T}.")
    print(f"The calculated integer value for D is {D}.")
    
    # Determine the best packing orientation for the final print statement
    if num_blocks == n_blocks_1:
        n_x, n_y, n_z = math.floor(C/T), math.floor(C/T), math.floor(C/D)
        dim_x, dim_y, dim_z = T, T, D
    elif num_blocks == n_blocks_2:
        n_x, n_y, n_z = math.floor(C/T), math.floor(C/D), math.floor(C/T)
        dim_x, dim_y, dim_z = T, D, T
    else: # num_blocks == n_blocks_3
        n_x, n_y, n_z = math.floor(C/D), math.floor(C/T), math.floor(C/T)
        dim_x, dim_y, dim_z = D, T, T

    final_calc_string = f"2 * (floor({C}/{dim_x}) * floor({C}/{dim_y}) * floor({C}/{dim_z})) = 2 * ({n_x} * {n_y} * {n_z}) = {num_prisms}"
    
    print("\nThe number of triangular prisms that can fit completely within the cube is calculated by packing pairs of prisms (blocks of size TxTxD).")
    print("The final calculation for the maximum number of prisms is:")
    print(final_calc_string)
    
    return num_prisms

# Execute the function and print the final answer in the required format
final_answer = solve_puzzle()
print(f"\n<<< {final_answer} >>>")