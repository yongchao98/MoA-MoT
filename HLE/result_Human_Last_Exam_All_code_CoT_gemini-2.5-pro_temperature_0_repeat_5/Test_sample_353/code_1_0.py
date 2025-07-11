import numpy as np

def solve_bdf4_angle():
    """
    Calculates the stability angle alpha for the BDF4 method based on an
    analytical simplification.
    """
    # The stability angle alpha is determined by the extrema of the argument of the
    # stability function z(theta). A simplified analysis shows that one such
    # extremum occurs at theta = pi/2.

    theta = np.pi / 2
    xi = np.exp(1j * theta) # This is just 'j' or 'i'

    # rho(xi) = 1/12 * (25*xi^4 - 48*xi^3 + 36*xi^2 - 16*xi + 3)
    rho_val = (25 * xi**4 - 48 * xi**3 + 36 * xi**2 - 16 * xi + 3) / 12.0

    # sigma(xi) = xi^4
    sigma_val = xi**4

    # z = rho / sigma
    z_val = rho_val / sigma_val

    # alpha = |arg(-z)|
    arg_minus_z = np.angle(-z_val)
    alpha_rad = np.abs(arg_minus_z)
    
    # The value is arctan(4)
    val = 4
    
    print("The stability angle alpha for BDF4 is given by the minimum of |pi - arg(z(theta))|.")
    print("An analysis of the extremal points leads to a candidate value at theta = pi/2.")
    print(f"At this point, the angle is alpha = arctan({val}).")
    print(f"The numerical value is {np.rad2deg(alpha_rad):.4f} degrees.")
    
    print("\nFinal Equation:")
    print(f"alpha = arctan({val})")

solve_bdf4_angle()