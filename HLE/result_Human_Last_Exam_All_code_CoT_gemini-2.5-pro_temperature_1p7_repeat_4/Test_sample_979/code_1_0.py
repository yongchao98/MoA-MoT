def print_solution():
    """
    This function prints the derived symbolic expressions for the magnetic field
    H inside and outside the spherical shell, and states the correct answer choice.
    """
    
    # Define the expressions for the H-field in both regions
    # These are derived from first principles using magnetostatics
    h_in_expression = "H_in(r, theta) = (2 * mu_0 / mu) * (K_0 / (1 + (2 * mu_0 / mu))) * z_hat"
    h_out_expression = "H_out(r, theta) = (K_0 / (1 + (2 * mu_0 / mu))) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"
    
    print("Based on the derivation using the magnetic scalar potential and boundary conditions:")
    print("-" * 70)

    # Print the field inside the sphere
    print("The magnetic field inside the sphere (0 < r < R) is uniform and given by:")
    print(f"  {h_in_expression}\n")
    print("where:")
    print("  K_0: current amplitude")
    print("  mu_0: permeability of free space")
    print("  mu: permeability of the material inside the sphere")
    print("  z_hat: unit vector in the z-direction")
    print("-" * 70)

    # Print the field outside the sphere
    print("The magnetic field outside the sphere (r > R) is a dipole field given by:")
    print(f"  {h_out_expression}\n")
    print("where:")
    print("  R: radius of the sphere")
    print("  r, theta: spherical coordinates")
    print("  r_hat, theta_hat: unit vectors in spherical coordinates")
    print("-" * 70)
    
    # State the conclusion
    print("These expressions for the magnetic field in both regions match Answer Choice E.")

# Execute the function to display the results
print_solution()