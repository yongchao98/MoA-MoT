def print_magnetic_field_solution():
    """
    Prints the mathematical expressions for the magnetic field H in both regions,
    based on the derived solution which corresponds to option B.
    """

    region_1_formula = (
        "In the region 0 < r < R_p:\n"
        "  H = M_0 * ((2*R_p**3 + R**3) / (3*R**3)) * "
        "(-cos(theta) * r_hat + sin(theta) * theta_hat)"
    )

    region_2_radial_component = (
        "- (2*M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * r_hat"
    )

    region_2_theta_component = (
        "+ (M_0 / 3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * theta_hat"
    )

    region_2_formula = (
        "In the region R_p < r < R:\n"
        f"  H = {region_2_radial_component} \n"
        f"      {region_2_theta_component}"
    )

    final_expression_H1 = "H(r < R_p) = M_0 * (2*R_p**3 + R**3)/(3*R**3) * [-cos(theta) i_r + sin(theta) i_theta]"
    final_expression_H2 = "H(R_p < r < R) = " \
                        "- (2*M_0/3) * [(R_p/R)**3 - (R_p/r)**3] * cos(theta) * i_r " \
                        "+ (M_0/3) * [2*(R_p/R)**3 + (R_p/r)**3] * sin(theta) * i_theta"


    print("The derived magnetic field H(r, theta) is:")
    print("\nIn the region 0 < r < R_p:")
    print("  H = M_0 * ( (2 * R_p^3 + R^3) / (3 * R^3) ) * ( -cos(theta) i_r + sin(theta) i_theta )")
    
    print("\nIn the region R_p < r < R:")
    print("  H = - (2*M_0/3) * [ (R_p/R)^3 - (R_p/r)^3 ] * cos(theta) i_r + (M_0/3) * [ 2*(R_p/R)^3 + (R_p/r)^3 ] * sin(theta) i_theta")

print_magnetic_field_solution()