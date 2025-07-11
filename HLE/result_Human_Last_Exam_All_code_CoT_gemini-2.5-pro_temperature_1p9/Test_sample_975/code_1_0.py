def solve_magnetostatic_shielding():
    """
    This function prints the derived expressions for the magnetic field H
    in the two specified regions.
    """

    # Expressions for the magnetic field H derived from the analytical solution.
    # The notation used is:
    # M_0: Magnitude of the uniform magnetization
    # R_p: Radius of the magnetized sphere
    # R: Radius of the outer conducting cavity
    # r, theta: Spherical coordinates
    # i_r, i_theta: Unit vectors in spherical coordinates

    h_field_inside = "H = M_0 * (2*R_p**3 + R**3) / (3*R**3) * (-cos(theta) * i_r + sin(theta) * i_theta)"
    
    h_field_outside_r = "H_r = - (2*M_0 / 3) * ( (R_p/R)**3 - (R_p/r)**3 ) * cos(theta) * i_r"
    h_field_outside_theta = "H_theta = (M_0 / 3) * ( 2*(R_p/R)**3 + (R_p/r)**3 ) * sin(theta) * i_theta"
    h_field_outside = f"{h_field_outside_r} + {h_field_outside_theta}"

    print("The derived magnetic field H is:")
    print("-" * 30)
    print("In the region 0 < r < R_p:")
    print(f"  {h_field_inside}")
    print("\nWhich corresponds to the expression:")
    print("  H = M_0 * ((2*R_p^3 + R^3) / (3*R^3)) * (-cos(theta) i_r + sin(theta) i_theta)")

    print("\nIn the region R_p < r < R:")
    # Combining the two components for a more readable full vector expression
    print("  H = [-2*M_0/3 * ((R_p/R)^3 - (R_p/r)^3) * cos(theta)] i_r")
    print("      + [M_0/3 * (2*(R_p/R)^3 + (R_p/r)^3) * sin(theta)] i_theta")

solve_magnetostatic_shielding()