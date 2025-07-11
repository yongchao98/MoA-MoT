def print_final_equation():
    """
    This function prints the derived equations for the magnetic field H
    in both regions, using symbolic variable names.
    """
    # Symbolic representation of variables
    K0 = "K_0"
    mu = "μ"
    mu0 = "μ_0"
    R = "R"
    r = "r"
    theta = "θ"
    z_hat = "ẑ"
    r_hat = "r̂"
    theta_hat = "θ̂"

    print("The final derived magnetic field H is:")
    print("H(r, θ) =")

    # Print H_in for r < R
    h_in_coeff_part1 = f"(2 * {mu0}) / {mu}"
    h_in_coeff_part2 = f"{K0} / (1 + (2 * {mu0}) / {mu})"
    print(f"  for 0 < r < R:  ({h_in_coeff_part1}) * ({h_in_coeff_part2}) * {z_hat}")

    # Print H_out for r > R
    h_out_coeff = f"{K0} / (1 + (2 * {mu0}) / {mu})"
    h_out_radial_part = f"({R}³ / {r}³)"
    h_out_angular_part = f"(2 * cos({theta}) * {r_hat} + sin({theta}) * {theta_hat})"
    print(f"  for R < r < ∞:  ({h_out_coeff}) * {h_out_radial_part} * {h_out_angular_part}")

    print("\nThis result matches answer choice E.")

print_final_equation()