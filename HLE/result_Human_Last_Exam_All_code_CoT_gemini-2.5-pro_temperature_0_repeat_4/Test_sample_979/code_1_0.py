def solve_and_print_magnetic_field():
    """
    This function prints the symbolic expressions for the magnetic field H
    inside and outside the described spherical shell.
    The expressions are derived from first principles of magnetostatics.
    """

    # Define symbolic constants as strings for printing
    K0 = "K_0"
    mu = "μ"
    mu0 = "μ_0"
    R = "R"
    r = "r"
    theta = "θ"

    # The derived expression for the H-field inside the sphere (0 < r < R)
    # This is a uniform field in the z-direction.
    # The expression is (2 * mu_0 / mu) * (K_0 / (1 + 2 * mu_0 / mu)) * z_hat
    # We print out each number in the equation as requested.
    H_in_numerator_coeff = 2
    H_in_denominator_const = 1
    H_in_denominator_coeff = 2
    
    H_in_str = (f"({H_in_numerator_coeff} * {mu0} / {mu}) * "
                f"({K0} / ({H_in_denominator_const} + {H_in_denominator_coeff} * {mu0} / {mu})) * z_hat")

    # The derived expression for the H-field outside the sphere (r > R)
    # This is a dipole field.
    # The expression is (K_0 / (1 + 2 * mu_0 / mu)) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    H_out_denominator_const = 1
    H_out_denominator_coeff = 2
    H_out_angular_coeff_r = 2
    
    H_out_str = (f"({K0} / ({H_out_denominator_const} + {H_out_denominator_coeff} * {mu0} / {mu})) * "
                 f"({R}^3 / {r}^3) * ({H_out_angular_coeff_r}*cos({theta})*r_hat + sin({theta})*theta_hat)")

    # Print the final results in a structured format
    print("The derived magnetic field H(r, θ) is:")
    print("-------------------------------------------------")
    print("For 0 < r < R (inside the sphere):")
    print(f"  H_in(r, θ) = {H_in_str}")
    print("\nFor r > R (outside the sphere):")
    print(f"  H_out(r, θ) = {H_out_str}")
    print("-------------------------------------------------")
    print("\nThese expressions correspond to Answer Choice E.")

solve_and_print_magnetic_field()