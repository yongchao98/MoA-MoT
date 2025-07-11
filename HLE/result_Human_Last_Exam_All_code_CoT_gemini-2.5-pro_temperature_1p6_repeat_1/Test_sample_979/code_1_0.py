def solve_and_print_magnetic_field():
    """
    This function prints the derived analytical expressions for the magnetic field
    H(r, theta) inside and outside the spherical shell.
    """

    # Define common variables and parts of the equations as strings
    common_denominator = "(1 + 2 * mu_0 / mu)"
    k0_str = "K_0"

    # --- Inside the sphere (0 < r < R) ---
    # Coefficient for the inside field
    h_in_coeff = f"(2 * mu_0 / mu) * {k0_str} / {common_denominator}"
    # The field is uniform in the z-direction
    h_in_field = f"{h_in_coeff} * z_hat"

    # --- Outside the sphere (r > R) ---
    # Coefficient for the outside field
    h_out_coeff_part1 = f"{k0_str} / {common_denominator}"
    h_out_coeff_part2 = "R^3 / r^3"
    h_out_coeff = f"{h_out_coeff_part1} * {h_out_coeff_part2}"
    # The field has a dipole form
    h_out_vector = "(2 * cos(theta) * r_hat + sin(theta) * theta_hat)"
    h_out_field = f"{h_out_coeff} * {h_out_vector}"
    
    # Print the final solution
    print("The derived magnetic field H(r, theta) is:")
    print("For 0 < r < R (inside the sphere):")
    print(f"H_in = {h_in_field}\n")
    
    print("For r > R (outside the sphere):")
    print(f"H_out = {h_out_field}\n")

    print("These expressions correspond to Answer Choice E.")

solve_and_print_magnetic_field()