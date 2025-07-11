def solve_magnetic_field():
    """
    This function prints the derived equations for the magnetic field H
    in the two specified regions, based on the analytical solution.
    """
    
    # --- Introduction to the results ---
    print("The magnetic field H(r, theta) is determined by solving Laplace's equation for the magnetic scalar potential")
    print("and applying the appropriate boundary conditions at r = R_p and r = R.")
    print("The final expressions for the magnetic field in each region are:")
    
    # --- Region 1: 0 < r < R_p ---
    print("\nIn the region 0 < r < R_p:")
    # The coefficient is M_0 * (2*R_p^3 + R^3) / (3*R^3)
    # The vector part is (-cos(theta) i_r + sin(theta) i_theta)
    # This matches Choice B.
    print("  H = M_0 * (2*R_p**3 + R**3) / (3*R**3) * (-cos(theta) i_r + sin(theta) i_theta)")

    # --- Region 2: R_p < r < R ---
    print("\nIn the region R_p < r < R:")
    # The H field is broken into its r and theta components as in the answer choices.
    # H_r = -2/3 * M_0 * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r
    # H_theta = 1/3 * M_0 * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta
    # This also matches Choice B.
    print("  H = - (2*M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r")
    print("      + (M_0 / 3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta")

# Execute the function to display the answer
solve_magnetic_field()