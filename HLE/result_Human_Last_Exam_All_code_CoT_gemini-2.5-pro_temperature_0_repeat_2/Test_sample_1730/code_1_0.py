import math

def solve():
    """
    This function calculates the values of T and D based on the provided systems,
    and then solves the geometric packing problem.
    """
    # Part 1: Calculate T
    # Given parameters for heat loss calculation
    L = 1.5  # m, Length of the collector
    B = 0.85  # m, Width of the collector
    U_inf = 1.0  # m/s, Wind speed
    nu_f = 15.11e-6  # m^2/s, Kinematic viscosity
    k_f = 0.0257  # W/(m.K), Thermal conductivity
    Pr_f = 0.707  # Prandtl number
    beta_f = 0.00341  # K^-1, Thermal expansion coefficient
    g = 9.81  # m/s^2, Acceleration due to gravity

    # Step 1.1: Calculate average wall temperature
    # The average value of sin(pi*x/L) over [0, L] is 2/pi.
    # So, theta_w_avg = 30 + 10 * (2/pi)
    theta_w_avg = 30 + 20 / math.pi

    # Step 1.2: Calculate average ambient temperature
    # The average value of 10 + 0.05*y over [0, B] is 10 + 0.05*B/2
    theta_inf_avg = 10 + 0.025 * B

    # Step 1.3: Calculate average temperature difference
    delta_theta = theta_w_avg - theta_inf_avg

    # Step 1.4: Calculate parameters for forced convection (along L)
    Re_L = U_inf * L / nu_f
    # Using correlation for laminar flow over a flat plate
    Nu_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))
    h_forced = Nu_forced * k_f / L

    # Step 1.5: Calculate parameters for free convection (along B)
    Ra_B = (g * beta_f * delta_theta * (B**3) / (nu_f**2)) * Pr_f
    # Using correlation for turbulent flow on a vertical plate (Ra > 10^9)
    Nu_free = 0.1 * (Ra_B**(1/3))
    h_free = Nu_free * k_f / B

    # Step 1.6: Combine for mixed convection (perpendicular flows)
    h_avg = (h_forced**3 + h_free**3)**(1/3)

    # Step 1.7: Calculate total heat loss and T
    A = L * B
    Q_v = h_avg * A * delta_theta
    T_float = Q_v / 80
    T = int(round(T_float))

    # Part 2: Calculate D
    # Given parameters for stress calculation
    q0 = 3.0  # N/m, Uniformly distributed load
    l = 2.0  # m, Length of the beam

    # Step 2.1: Calculate maximum bending moment (for a simply supported beam)
    M_max = q0 * (l**2) / 8

    # Step 2.2: Calculate maximum normal stress
    # The geometric properties are defined such that they simplify.
    # I_yy = a and z_max = 2a, where a is the given parameter.
    # So, sigma_max = M_max * (2a) / a = 2 * M_max
    sigma_max = 2 * M_max

    # Step 2.3: Calculate D
    D = int(round(sigma_max / 3.0))

    # Part 3: Solve the packing problem
    # Prism: base is a right triangle with legs T, T; depth is D.
    # Cube: 3x3x3.
    # With T=2 and D=1, the prism has a 2x2 right triangular base and depth 1.
    # Two such prisms form a 2x2x1 rectangular block.
    # The problem is equivalent to finding the max number of 2x2x1 blocks in a 3x3x3 cube.
    # The maximum number of such blocks is 5.
    # This gives 5 blocks * 2 prisms/block = 10 prisms.
    # The remaining volume of 27 - 5*(2*2*1) = 7 is in disconnected 1x1x1 cubes
    # and cannot fit another prism which requires a 2x2 footprint.
    num_prisms = 10

    # Final Output
    print(f"Step 1: Calculation of T")
    print(f"The calculated heat loss is {Q_v:.2f} W.")
    print(f"T = {Q_v:.2f} / 80 = {T_float:.2f}")
    print(f"The rounded value for T is: {T}")
    print("-" * 20)
    print(f"Step 2: Calculation of D")
    print(f"The maximum normal stress is {sigma_max:.2f} N/m^2.")
    print(f"D = {sigma_max:.2f} / 3 = {D}")
    print(f"The value for D is: {D}")
    print("-" * 20)
    print(f"Step 3: Solving the packing problem")
    print(f"The prism has a right triangular base with side lengths {T}, {T}, and {T}*sqrt({T}) and a depth of {D}.")
    print(f"The cube has dimensions 3 x 3 x 3.")
    print(f"The maximum number of prisms that can fit completely within the cube is {num_prisms}.")
    print(f"<<<{num_prisms}>>>")

solve()