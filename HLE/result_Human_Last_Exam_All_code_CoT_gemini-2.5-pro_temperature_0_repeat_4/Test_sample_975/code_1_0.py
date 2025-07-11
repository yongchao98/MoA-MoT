def print_solution():
    """
    Prints the final equations for the magnetic field H in both regions.
    """
    
    # The solution corresponds to option B
    
    # Expression for the magnetic field inside the shield (0 < r < R_p)
    h_field_inside = r"In the region 0 < r < R_p:"
    h_eq_inside = r"  H = M_0 * ((2*R_p**3 + R**3) / (3*R**3)) * (-cos(theta) * i_r + sin(theta) * i_theta)"
    
    # Expression for the magnetic field between the shield and the conductor (R_p < r < R)
    h_field_outside = r"In the region R_p < r < R:"
    h_eq_outside_r = r"  H_r = -(2*M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r"
    h_eq_outside_theta = r"  H_theta = (M_0 / 3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta"
    h_eq_outside_full = r"  H = -" + r"2*M_0" + r"/3 * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r + " + r"M_0" + r"/3 * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta"


    print(h_field_inside)
    # To present the equation clearly, we print each term of the final equation
    print("  H = M_0 * ( (2*R_p^3 + R^3) / (3*R^3) ) * ( -cos(theta) i_r + sin(theta) i_theta )")
    print("\n" + h_field_outside)
    print("  H = -(2*M_0/3) * [ (R_p/R)^3 - (R_p/r)^3 ] * cos(theta) i_r + (M_0/3) * [ 2*(R_p/R)^3 + (R_p/r)^3 ] * sin(theta) i_theta")
    
    # To be more explicit as requested by the prompt, let's print each number in the final equation
    print("\nDecomposition of the final equations:")
    print("Region 1 (0 < r < R_p):")
    print("  Coefficient: M_0 * ( (2*R_p^3 + R^3) / (3*R^3) )")
    print("  Vector part: ( -cos(theta) i_r + sin(theta) i_theta )")
    print("  Numbers in equation: 2, 3, 3, 3")
    
    print("\nRegion 2 (R_p < r < R):")
    print("  Radial component (H_r):")
    print("    Coefficient: -(2*M_0/3)")
    print("    Spatial part: [ (R_p/R)^3 - (R_p/r)^3 ]")
    print("    Angular part: cos(theta) i_r")
    print("    Numbers in equation: 2, 3, 3, 3")
    print("  Theta component (H_theta):")
    print("    Coefficient: (M_0/3)")
    print("    Spatial part: [ 2*(R_p/R)^3 + (R_p/r)^3 ]")
    print("    Angular part: sin(theta) i_theta")
    print("    Numbers in equation: 3, 2, 3, 3")


print_solution()