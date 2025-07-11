def print_magnetic_field_solution():
    """
    This function prints the final derived expressions for the magnetic field H
    in the two specified regions of the substation shielding problem.
    """

    print("The magnetic field H(r, theta) is determined in two regions:")
    print("=" * 65)

    # Region 1: Inside the shield (0 < r < Rp)
    # H_1 = M_0 * ( (R**3 + 2*R_p**3)/(3*R**3) ) * (-cos(theta) i_r + sin(theta) i_theta)
    # The numerical coefficients involved are 2, 3, -1, 1.
    print("In the region 0 < r < R_p:")
    print("  H = M_0 * (R**3 + 2*R_p**3)/(3*R**3) * (-1*cos(theta) i_r + 1*sin(theta) i_theta)")
    print("")

    # Region 2: Between the shield and the conductor (Rp < r < R)
    # H_2 = - (2*M_0/3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) i_r
    #         + (M_0/3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) i_theta
    # The numerical coefficients are -2, 3, 2, 1.
    print("In the region R_p < r < R:")
    radial_component = "- (2*M_0/3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) i_r"
    angular_component = "+ (1*M_0/3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) i_theta"
    print(f"  H = {radial_component}")
    print(f"      {angular_component}")
    
    print("=" * 65)

print_magnetic_field_solution()