def print_magnetic_field_solution():
    """
    This function prints the derived magnetic field H(r, theta)
    for the given spherical shell problem.
    """
    
    # The solution is divided into two regions: inside the sphere (r < R)
    # and outside the sphere (r > R).

    # Expression for the magnetic field inside the sphere (0 < r < R)
    # The field is uniform and directed along the z-axis.
    # H_in = (2 * mu_0 / mu) * (K_0 / (1 + (2 * mu_0 / mu))) * z_hat
    inside_field = "H(r, theta) = (2 * mu_0 / mu) * (K_0 / (1 + (2 * mu_0 / mu))) * z_hat"

    # Expression for the magnetic field outside the sphere (r > R)
    # The field has the form of a magnetic dipole.
    # H_out = (K_0 / (1 + (2 * mu_0 / mu))) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)
    outside_field = "H(r, theta) = (K_0 / (1 + (2 * mu_0 / mu))) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)"

    print("The derived magnetic field H(r, theta) is:")
    print("-" * 40)
    print("For the region 0 < r < R (inside the sphere):")
    print(inside_field)
    print("\nFor the region r > R (outside the sphere):")
    print(outside_field)
    print("-" * 40)
    print("\nThis corresponds to Answer Choice E.")

# Execute the function to print the solution
print_magnetic_field_solution()