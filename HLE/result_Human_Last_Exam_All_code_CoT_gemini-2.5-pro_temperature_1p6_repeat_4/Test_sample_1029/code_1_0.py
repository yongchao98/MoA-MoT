def display_poynting_vector():
    """
    This script computes and displays the formula for the Poynting vector S
    for a moving, charged cylindrical rod in an external electric field.
    The result is presented symbolically.
    """

    print("The Poynting vector S is calculated for two regions: inside and outside the rod.")
    print("A cylindrical coordinate system (r, phi, z) is used.")
    print("The unit vectors are r_hat (radial), phi_hat (azimuthal), and k_hat (axial).\n")

    # The equation has two main components: one in the axial direction (k_hat)
    # and one in the radial direction (r_hat). The numbers in the equation
    # are the numerical coefficients 2 and 4 in the denominators.

    # Inside the rod (r <= R)
    print("--- Inside the rod (r <= R) ---")
    s_in_str = "S_in = (rho**2 * v * r**2 / (4 * epsilon_0)) * k_hat - (E * rho * v * r / 2) * r_hat"
    print(s_in_str)
    print("In this equation, the coefficients are:")
    print("  k_hat component: 1/4")
    print("  r_hat component: -1/2\n")

    # Outside the rod (r > R)
    print("--- Outside the rod (r > R) ---")
    s_out_str = "S_out = (rho**2 * v * R**4 / (4 * epsilon_0 * r**2)) * k_hat - (E * rho * v * R**2 / (2 * r)) * r_hat"
    print(s_out_str)
    print("In this equation, the coefficients are:")
    print("  k_hat component: 1/4")
    print("  r_hat component: -1/2\n")


    print("--- Variable Definitions ---")
    print("  rho       : uniform volume charge density")
    print("  E         : magnitude of the external axial electric field")
    print("  v         : speed of the rod along its axis")
    print("  R         : radius of the rod")
    print("  r         : radial distance from the rod's axis")
    print("  epsilon_0 : permittivity of free space")


if __name__ == "__main__":
    display_poynting_vector()
