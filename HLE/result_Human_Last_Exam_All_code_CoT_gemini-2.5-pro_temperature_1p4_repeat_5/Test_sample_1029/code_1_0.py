def solve_poynting_vector():
    """
    This script calculates and prints the symbolic expression for the Poynting vector
    both inside and outside the described cylindrical rod.
    
    The final equations for the vector components will be printed.
    """
    
    # --- Inside the rod (for r <= R) ---
    print("Poynting vector inside the rod (for a radial distance r <= R):")
    print("In cylindrical coordinates, the Poynting vector S has a radial (r) and an axial (z) component.")
    
    # Equation for the radial component inside
    # S_r = - (E * rho * v * r) / 2
    print("\nThe equation for the radial component (S_r) is:")
    print("S_r = - (E * rho * v * r) / 2")

    # Equation for the axial component inside
    # S_z = (rho^2 * v * r^2) / (4 * epsilon_0)
    print("\nThe equation for the axial component (S_z) is:")
    print("S_z = (rho^2 * v * r^2) / (4 * epsilon_0)")

    print("\nThus, the final vector equation inside the rod is:")
    print("S_vector = [ - (E * rho * v * r) / 2 ] r_hat + [ (rho^2 * v * r^2) / (4 * epsilon_0) ] z_hat")

    print("\n" + "="*70 + "\n")

    # --- Outside the rod (for r > R) ---
    print("Poynting vector outside the rod (for a radial distance r > R):")

    # Equation for the radial component outside
    # S_r = - (E * rho * v * R^2) / (2 * r)
    print("\nThe equation for the radial component (S_r) is:")
    print("S_r = - (E * rho * v * R^2) / (2 * r)")

    # Equation for the axial component outside
    # S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)
    print("\nThe equation for the axial component (S_z) is:")
    print("S_z = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)")
    
    print("\nThus, the final vector equation outside the rod is:")
    print("S_vector = [ - (E * rho * v * R^2) / (2 * r) ] r_hat + [ (rho^2 * v * R^4) / (4 * epsilon_0 * r^2) ] z_hat")


# Execute the function to print the solution
solve_poynting_vector()
