def solve_magnetic_field():
    """
    This function defines and prints the symbolic expressions for the magnetic field
    H inside and outside the spherical shell, based on the analytical solution.
    """

    # --- Symbolic variable names ---
    K0 = "K_0"
    mu = "mu"
    mu0 = "mu_0"
    R = "R"
    r = "r"
    theta = "theta"
    r_hat = "r_hat"
    theta_hat = "theta_hat"
    z_hat = "z_hat"

    # --- Construct the expression for H-field inside the sphere (0 < r < R) ---
    # This corresponds to the term in option E.
    # H_in = (2 * mu_0 / mu) * (K_0 / (1 + 2 * mu_0 / mu)) * z_hat
    H_in_coeff_part1 = f"(2 * {mu0} / {mu})"
    H_in_coeff_part2 = f"({K0} / (1 + 2 * {mu0} / {mu}))"
    H_in_expression = f"{H_in_coeff_part1} * {H_in_coeff_part2} * {z_hat}"

    # --- Construct the expression for H-field outside the sphere (r > R) ---
    # This corresponds to the term in option E.
    # H_out = (K_0 / (1 + 2 * mu_0 / mu)) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    H_out_coeff_part1 = f"({K0} / (1 + 2 * {mu0} / {mu}))"
    H_out_spatial_dependence = f"({R}^3 / {r}^3)"
    H_out_vector_part = f"(2 * cos({theta}) * {r_hat} + sin({theta}) * {theta_hat})"
    H_out_expression = f"{H_out_coeff_part1} * {H_out_spatial_dependence} * {H_out_vector_part}"

    # --- Print the final results ---
    print("The magnetic field H(r, theta) is given by:")
    print("\nFor 0 < r < R (inside the sphere):")
    print(f"  H_in = {H_in_expression}")
    print("\nFor r > R (outside the sphere):")
    print(f"  H_out = {H_out_expression}")

# Execute the function to print the solution
solve_magnetic_field()