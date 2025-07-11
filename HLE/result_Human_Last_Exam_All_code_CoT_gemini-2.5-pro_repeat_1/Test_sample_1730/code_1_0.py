import numpy as np
from scipy.integrate import quad

def solve_and_print():
    """
    This function calculates the values of T and D based on the problem descriptions,
    and then determines how many triangular prisms can fit into a cube.
    """
    # Part 1: Calculation of T from the heat transfer problem
    # Given constants
    L = 1.5         # m
    B = 0.85        # m
    U_inf = 1.0     # m/s
    k_f = 0.0257    # W/(m.K)
    nu_f = 15.11e-6 # m^2/s
    Pr_f = 0.707    # Prandtl number

    # The local heat transfer coefficient h(x) is h(x) = C_h * x^(-1/2)
    # Calculate the constant C_h
    C_h = 0.332 * k_f * (Pr_f**(1/3)) * ((U_inf / nu_f)**0.5)

    # The total heat loss Q is the integral of h(x) * (theta_w(x) - theta_inf(y)) over the area.
    # We first integrate (theta_w(x) - theta_inf(y)) with respect to y from 0 to B.
    # integral_y = integral_0^B [ (30 + 10sin(pi*x/L)) - (10 + 0.05y) ] dy
    #            = [ (20 + 10sin(pi*x/L))y - 0.025y^2 ] from 0 to B
    #            = (20 + 10sin(pi*x/L))B - 0.025B^2
    # This result is a function of x, let's call it G(x).
    # The total heat loss Q is the integral of h(x)*G(x) from 0 to L.
    def integrand_for_Q(x, B, L, C_h):
        # This is the function h(x) * G(x)
        if x == 0:
            return 0.0 # The integral is convergent, but the function diverges at x=0
        G_x = (20 * B) + (10 * B * np.sin(np.pi * x / L)) - (0.025 * B**2)
        h_x = C_h * x**(-0.5)
        return h_x * G_x

    # Perform numerical integration to find the total heat loss Q_dot_V
    Q_dot_V, _ = quad(integrand_for_Q, 0, L, args=(B, L, C_h))

    # Calculate T by dividing Q_dot_V by 80 and rounding to the nearest integer
    T_val = int(round(Q_dot_V / 80.0))

    # Part 2: Calculation of D from the beam bending problem
    # Given constants
    q0 = 3.0  # N/m
    l = 2.0   # m

    # For a simply supported beam with a uniform load, max bending moment is at the center
    M_max = (q0 * l**2) / 8.0

    # The problem is set up such that I_y simplifies to 'a' and z_max is '2a'.
    # sigma_xx_max = (M_max * z_max) / I_y = (M_max * 2a) / a = 2 * M_max
    sigma_xx_max = 2 * M_max

    # Calculate D by dividing sigma_xx_max by 3
    D_val = int(round(sigma_xx_max / 3.0))

    # Part 3: Solving the packing problem
    # The problem statement for the prism's base (right triangle with sides T, T, T*sqrt(T))
    # is contradictory for T=1. We assume the intended base is a right isosceles triangle
    # with legs of length T.
    
    # With this assumption, two prisms can be joined to form a T x T x D cuboid.
    # We can calculate how many such cuboids fit in the 3x3x3 cube.
    cube_dim = 3
    prisms_per_cuboid = 2

    # Calculate the total number of prisms
    num_prisms = (cube_dim / T_val) * (cube_dim / T_val) * (cube_dim / D_val) * prisms_per_cuboid

    # Print the final answer in the requested equation format
    print("The final calculation for the number of prisms is:")
    print(f"({cube_dim} / {T_val}) * ({cube_dim} / {T_val}) * ({cube_dim} / {D_val}) * {prisms_per_cuboid} = {int(num_prisms)}")
    
solve_and_print()
<<<54>>>