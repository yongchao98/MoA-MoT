def display_magnetic_field_equations():
    """
    This function prints the equations for the magnetic field H in the two regions,
    corresponding to the correct answer choice.
    """
    
    print("The magnetic field H(r, theta) is determined in two regions:")
    print("-" * 60)

    # Display the equation for the region inside the shield (0 < r < R_p)
    print("In the region 0 < r < R_p:")
    # H = M0 * (2*Rp^3 + R^3)/(3*R^3) * (-cos(theta) i_r + sin(theta) i_theta)
    print("  H = M_0 * ((2*R_p**3 + R**3) / (3*R**3)) * (-cos(theta) i_r + sin(theta) i_theta)")
    print("\nWhich can be written as:")
    print("  H_r = -M_0 * ((2*R_p**3 + R**3) / (3*R**3)) * cos(theta)")
    print("  H_theta = M_0 * ((2*R_p**3 + R**3) / (3*R**3)) * sin(theta)")
    
    print("\n" + "-" * 60 + "\n")

    # Display the equation for the region between the shield and the conductor (R_p < r < R)
    print("In the region R_p < r < R:")
    # H = -2*M0/3 * [ (Rp/R)^3 - (Rp/r)^3 ] * cos(theta) i_r + M0/3 * [ 2*(Rp/R)^3 + (Rp/r)^3 ] * sin(theta) i_theta
    print("  H_r = - (2*M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta)")
    print("  H_theta = (M_0 / 3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta)")

display_magnetic_field_equations()