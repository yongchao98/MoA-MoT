def display_electric_field_solution():
    """
    This function prints the symbolic expressions for the electric field
    in the specified regions, as derived from solving Laplace's equation
    with the given boundary conditions.
    """
    # Symbolic variable names
    P0 = "P_0"
    eps0 = "varepsilon_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    theta = "theta"
    r_hat = "r_hat"
    theta_hat = "theta_hat"

    # --- Electric Field inside the sensor (r < R_p) ---
    print("The electric field for r < R_p (inside the sensor) is:")
    
    # Breaking down the expression for clarity
    # E_in = - (Coefficient) * (Vector Part)
    # The term (cos(theta)*r_hat - sin(theta)*theta_hat) is the unit vector z_hat.
    # The field is uniform and directed opposite to the polarization.
    
    coeff_in = f"({P0} / (3 * {eps0}))"
    factor_in = f"(1 - ({Rp}/{R})^3)"
    vector_part_z = f"(cos({theta}) * {r_hat} - sin({theta}) * {theta_hat})"
    
    print(f"E_in = - [{coeff_in}] * [{factor_in}] * [{vector_part_z}]\n")
    
    # --- Electric Field in the free space (R_p < r < R) ---
    print("The electric field for R_p < r < R (in the free space) is:")
    
    # This field is a superposition of two fields:
    # 1. A uniform field due to the charges induced on the outer conducting shell.
    # 2. A dipole field due to the polarized sphere.
    
    # Term 1: Uniform field
    coeff_out1 = f"({P0} / (3 * {eps0}))"
    factor_out1 = f"({Rp}/{R})^3"
    vector_part1 = f"(cos({theta}) * {r_hat} - sin({theta}) * {theta_hat})"
    term1 = f"[{coeff_out1}] * [{factor_out1}] * [{vector_part1}]"
    
    # Term 2: Dipole field
    coeff_out2 = f"({P0} * {Rp}^3 / (3 * {eps0} * {r}^3))"
    vector_part2 = f"(2*cos({theta}) * {r_hat} + sin({theta}) * {theta_hat})"
    term2 = f"[{coeff_out2}] * [{vector_part2}]"
    
    print(f"E_out = {term1} + {term2}\n")

if __name__ == "__main__":
    display_electric_field_solution()
