import math

def solve():
    """
    This function solves the entire problem step-by-step.
    1. Calculates T from the heat transfer system.
    2. Calculates D from the beam mechanics system.
    3. Calculates how many prisms fit into the cube.
    """

    # Part 1: Calculation of T
    # Given constants for heat transfer
    L = 1.5  # m
    B = 0.85  # m
    U_inf = 1.0  # m/s
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m.K)
    Pr_f = 0.707
    g = 9.81  # m/s^2
    beta_f = 0.00341  # K^-1

    # Calculate average temperatures
    theta_w_avg = 30.0 + (2.0 * 10.0 / math.pi)
    theta_inf_avg = 10.0 + (0.05 * B / 2.0)
    delta_theta_avg = theta_w_avg - theta_inf_avg

    # Calculate dimensionless numbers for mixed convection
    Re_L = U_inf * L / nu_f
    Gr_L = (g * beta_f * delta_theta_avg * L**3) / (nu_f**2)
    Ra_L = Gr_L * Pr_f

    # Calculate Nusselt numbers for forced and natural convection
    # Flow is laminar for forced convection (Re_L < 5e5)
    Nu_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    # Flow is turbulent for natural convection (Ra_L > 1e9)
    Nu_natural = 0.1 * (Ra_L**(1/3))

    # Combine for mixed convection (perpendicular flows)
    Nu_mixed = (Nu_forced**3 + Nu_natural**3)**(1/3)

    # Calculate total heat loss Q_v
    h_avg = Nu_mixed * k_f / L
    A = L * B
    Q_v = h_avg * A * delta_theta_avg

    # Calculate T
    T_float = Q_v / 80.0
    T = int(round(T_float))

    # Part 2: Calculation of D
    # Given constants for beam mechanics
    q0 = 3.0  # N/m
    l = 2.0  # m

    # Calculate maximum bending moment for a simply supported beam with UDL
    M_max = q0 * l**2 / 8.0

    # The cross-section parameters are designed to simplify the stress calculation.
    # The moment of inertia I_yy simplifies to be numerically equal to 'a',
    # and Z_max is 2a. So, sigma_max = M_max * 2a / a = 2 * M_max.
    sigma_xx_max = 2.0 * M_max

    # Calculate D
    D_float = sigma_xx_max / 3.0
    D = int(round(D_float))

    # Part 3: Calculate how many prisms fit in the cube
    S = 3  # Cube side length

    # A right triangular prism with leg T and depth D can be paired with another
    # to form a rectangular block of dimensions T x T x D.
    # With T=2, the base sides are 2, 2, and 2*sqrt(2), which is a valid right triangle.
    block_dims = sorted([T, T, D], reverse=True)
    cube_dim = S

    # Calculate the number of blocks that can be packed
    num_blocks = math.floor(cube_dim / block_dims[0]) * \
                 math.floor(cube_dim / block_dims[1]) * \
                 math.floor(cube_dim / block_dims[2])
    
    # Total number of prisms is twice the number of blocks
    num_prisms = 2 * num_blocks
    
    # Print the final calculation as requested
    print(f"The calculated dimensions for the prism are T = {T} and D = {D}.")
    print("Two prisms are combined to form a rectangular block of size T x T x D.")
    print(f"The block dimensions are {T} x {T} x {D}.")
    print(f"The cube has a side length S = {S}.")
    print("\nThe number of prisms is found by packing these blocks into the cube:")
    print(f"Number of blocks = floor(S/T) * floor(S/T) * floor(S/D)")
    n_x = math.floor(S/T)
    n_y = math.floor(S/T)
    n_z = math.floor(S/D)
    print(f"Number of blocks = floor({S}/{T}) * floor({S}/{T}) * floor({S}/{D}) = {n_x} * {n_y} * {n_z} = {num_blocks}")
    print(f"Total number of prisms = 2 * (Number of blocks)")
    print(f"Final equation: {num_prisms} = 2 * {num_blocks}")

solve()
<<<6>>>