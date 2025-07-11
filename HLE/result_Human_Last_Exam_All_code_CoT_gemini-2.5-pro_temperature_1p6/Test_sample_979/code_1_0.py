def print_magnetic_field_solution():
    """
    This function prints the derived magnetic field H inside and outside the sphere.
    """
    
    # Define the components of the solution symbolically
    H_in_str = "H_in(r, theta) = (2 * mu_0 * K_0) / (mu * (1 + (2 * mu_0) / mu)) * z_hat"
    H_out_str = "H_out(r, theta) = K_0 / (1 + (2 * mu_0) / mu) * (R^3 / r^3) * (2 * cos(theta) * r_hat + sin(theta) * theta_hat)"

    print("The magnetic field H(r, theta) is found by solving Laplace's equation for the magnetic scalar potential with the given boundary conditions.")
    print("\nFor the region inside the sphere (0 < r < R):")
    print(H_in_str)
    
    print("\nFor the region outside the sphere (R < r < infinity):")
    print(H_out_str)
    
    print("\nThese expressions correspond to choice E.")
    print("\nHere are the final equations with all terms displayed explicitly:")
    
    # Breaking down the expression for clarity
    print("\n--- Inside the Sphere (0 < r < R) ---")
    print("H_in = [ 2 * (mu_0/mu) * K_0 / (1 + 2*(mu_0/mu)) ] * z_hat")
    
    print("\n--- Outside the Sphere (r > R) ---")
    print("H_out = [ K_0 / (1 + 2*(mu_0/mu)) ] * (R^3/r^3) * [ 2*cos(theta)*r_hat + sin(theta)*theta_hat ]")


# Execute the function to print the solution
print_magnetic_field_solution()