def print_magnetic_field_solution():
    """
    This function prints the derived expressions for the magnetic field H
    in the two specified regions, corresponding to the correct answer.
    The numbers in the equations are explicitly shown as per the instructions.
    """

    print("The magnetic field H(r, theta) is determined as follows:")

    # Solution for the region inside the shield (0 < r < R_p)
    print("\nIn the region 0 < r < R_p:")
    # The full vector expression is constructed from its components
    # The term ( - cos(theta) i_r + sin(theta) i_theta ) is a representation of the -i_z unit vector
    h1_coeff_numerator = "2*R_p**3 + R**3"
    h1_coeff_denominator = "3*R**3"
    vector_part_1 = "(-cos(theta) i_r + sin(theta) i_theta)"
    print(f"  H = M_0 * ({h1_coeff_numerator}) / ({h1_coeff_denominator}) * {vector_part_1}")


    # Solution for the region between the shield and the conductor (R_p < r < R)
    print("\nIn the region R_p < r < R:")
    # We print the radial and angular components separately for clarity
    hr2_coeff = f"-({2}*M_0/{3})"
    hr2_term = f"[(R_p/R)**{3} - (R_p/r)**{3}]"
    print(f"  H_r = {hr2_coeff} * {hr2_term} * cos(theta)")

    htheta2_coeff = f"(M_0/{3})"
    htheta2_term = f"[{2}*(R_p/R)**{3} + (R_p/r)**{3}]"
    print(f"  H_theta = {htheta2_coeff} * {htheta2_term} * sin(theta)")

# Execute the function to display the solution
print_magnetic_field_solution()