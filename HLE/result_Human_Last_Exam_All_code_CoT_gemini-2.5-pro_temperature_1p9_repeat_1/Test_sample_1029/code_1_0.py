def display_poynting_vector_formulas():
    """
    This function prints the derived symbolic formulas for the Poynting vector S
    for a moving, charged cylindrical rod in an external electric field.

    The vector S is presented in cylindrical coordinates (r, phi, z) with unit
    vectors r_hat, phi_hat, and k_hat.

    Variables used in the formulas:
    R: radius of the cylindrical rod
    rho: uniform volume charge density
    E: uniform external electric field along the axis (z-direction)
    v: speed of the rod along its axis (z-direction)
    r: radial distance from the center of the axis
    epsilon_0: permittivity of free space
    """
    print("The Poynting vector S describes the energy flux. Its formula depends on the location relative to the rod.")
    print("-" * 70)

    # Case 1: Inside the rod
    print("Case 1: Inside the rod (where radial distance r < R)")
    print("The Poynting vector S is the sum of a radial and an axial component: S = S_r * r_hat + S_z * k_hat\n")

    print("Radial Component (S_r):")
    print("S_r = -(E * rho * v * r) / 2")
    print("The final equation for the radial component consists of the following numbers and variables:")
    print("  - The coefficient: -1/2")
    print("  - The external electric field: E")
    print("  - The charge density: rho")
    print("  - The velocity: v")
    print("  - The radial distance: r\n")

    print("Axial Component (S_z):")
    print("S_z = (rho**2 * v * r**2) / (4 * epsilon_0)")
    print("The final equation for the axial component consists of the following numbers and variables:")
    print("  - The coefficient: 1/4")
    print("  - The charge density squared: rho^2")
    print("  - The velocity: v")
    print("  - The radial distance squared: r^2")
    print("  - This is all divided by the constant: epsilon_0")

    print("\n" + "=" * 70 + "\n")

    # Case 2: Outside the rod
    print("Case 2: Outside the rod (where radial distance r > R)")
    print("The Poynting vector S is also the sum of a radial and an axial component: S = S_r * r_hat + S_z * k_hat\n")

    print("Radial Component (S_r):")
    print("S_r = -(E * rho * v * R**2) / (2 * r)")
    print("The final equation for the radial component consists of the following numbers and variables:")
    print("  - The coefficient: -1/2")
    print("  - The external electric field: E")
    print("  - The charge density: rho")
    print("  - The velocity: v")
    print("  - The rod radius squared: R^2")
    print("  - This is all divided by the radial distance: r\n")

    print("Axial Component (S_z):")
    print("S_z = (rho**2 * v * R**4) / (4 * epsilon_0 * r**2)")
    print("The final equation for the axial component consists of the following numbers and variables:")
    print("  - The coefficient: 1/4")
    print("  - The charge density squared: rho^2")
    print("  - The velocity: v")
    print("  - The rod radius to the fourth power: R^4")
    print("  - This is all divided by the constant epsilon_0 and the radial distance squared r^2")
    
    print("-" * 70)

# Execute the function to print the solution
display_poynting_vector_formulas()