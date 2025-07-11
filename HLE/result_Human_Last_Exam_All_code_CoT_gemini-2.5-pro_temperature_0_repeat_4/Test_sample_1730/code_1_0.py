import numpy as np
from scipy.integrate import quad

def solve_puzzle():
    """
    This function calculates the values of T and D based on the physics problems,
    and then solves the geometric packing problem.
    """
    # Part 1: Calculate T
    # Constants for heat transfer
    L = 1.5  # m
    B = 0.85  # m
    U_inf = 1.0  # m/s
    nu_f = 15.11e-6  # m^2/s
    k_f = 0.0257  # W/(m.K)
    Pr_f = 0.707

    # The local heat transfer coefficient is h_x = C_h * x**(-0.5)
    # where C_h = 0.332 * k_f * Pr_f**(1/3) * (U_inf / nu_f)**(1/2)
    C_h = 0.332 * k_f * (Pr_f)**(1/3) * (U_inf / nu_f)**(0.5)

    # To find the total heat loss Q_V, we integrate the local heat flux over the plate area.
    # Q_V = integral from 0 to B of [ integral from 0 to L of h_x * (theta_w(x) - theta_inf(y)) dx ] dy
    # This can be separated and solved.

    # Integral of x**(-0.5) dx from 0 to L
    integral_x_term1 = 2 * L**0.5

    # Numerically integrate the sine term: integral of x**(-0.5) * sin(pi*x/L) dx from 0 to L
    def integrand_sin(x):
        # The function is integrable at x=0, but to avoid a 0/0 error in implementation,
        # we handle it as a special case. The limit of the integrand as x->0 is 0.
        if x == 0:
            return 0
        return x**(-0.5) * np.sin(np.pi * x / L)
    integral_x_term2, _ = quad(integrand_sin, 0, L)

    # Integral of (20 - 0.05*y) dy from 0 to B
    integral_y_term1 = 20 * B - 0.025 * B**2

    # Integral of 1 dy from 0 to B
    integral_y_term2 = B

    # Total heat loss Q_V is the sum of the integrated parts
    # Q_V = C_h * integral_x_term1 * integral_y_term1 + 10 * C_h * integral_x_term2 * integral_y_term2
    Q_V = C_h * (integral_x_term1 * integral_y_term1 + 10 * integral_x_term2 * integral_y_term2)

    # Calculate T
    T_float = Q_V / 80.0
    T = int(round(T_float))

    # Part 2: Calculate D
    # Constants for beam bending
    q0 = 3.0  # N/m
    l = 2.0  # m

    # Maximum bending moment for a simply supported beam with uniform load
    M_y_max = q0 * l**2 / 8.0

    # The moment of inertia I_yy and max distance Z_max are defined in terms of 'a'.
    # I_yy = a
    # Z_max = 2a
    # The 'a' term cancels out in the stress calculation.
    # sigma_xx_max = M_y_max * Z_max / I_yy = M_y_max * (2a) / a
    sigma_xx_max = 2 * M_y_max

    # Calculate D
    D = sigma_xx_max / 3.0

    # Part 3: Solve the packing problem
    # Cube dimensions
    cube_side = 3

    # Prism dimensions based on calculated T and D.
    # We assume the base is an isosceles right triangle with legs of length T.
    prism_leg = T
    prism_depth = int(D)

    # The volume of one prism
    prism_volume = 0.5 * prism_leg * prism_leg * prism_depth
    
    # The volume of the large cube
    cube_volume = cube_side**3

    # The number of prisms is the ratio of volumes, as they can be packed perfectly.
    num_prisms = int(cube_volume / prism_volume)

    # Print the final equation with the calculated numbers
    print(f"The number of prisms is calculated by dividing the cube's volume by the prism's volume.")
    print(f"Number of prisms = ({cube_side} * {cube_side} * {cube_side}) / (0.5 * {prism_leg} * {prism_leg} * {prism_depth}) = {num_prisms}")
    print(f"<<<{num_prisms}>>>")

solve_puzzle()