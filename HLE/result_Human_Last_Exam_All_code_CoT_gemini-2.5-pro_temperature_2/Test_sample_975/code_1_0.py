def display_magnetic_field_solution():
    """
    This function displays the symbolic solution for the magnetic field H in a two-region
    magnetostatic problem involving a magnetized sphere inside a conducting shell.
    The output corresponds to the correct answer choice.
    """

    print("The final expressions for the magnetic field H in the two regions are:")
    print("=" * 65)

    # --- Region 1: Inside the shield (0 < r < R_p) ---
    print("In the region 0 < r < R_p:")
    # H = M_0 * (2*R_p^3 + R^3) / (3*R^3) * (-cos(theta) i_r + sin(theta) i_theta)
    
    # We output each "number" (coefficient/term) in the equation.
    coeff_part = "M_0 * (2*R_p^3 + R^3) / (3*R^3)"
    vector_part = "(-cos(theta) * i_r + sin(theta) * i_theta)"
    
    print(f"  H = ({coeff_part}) * ({vector_part})")
    
    print("\n   Key numerical constants in the coefficient are 2 and 3.")
    
    print("\n" + "=" * 65 + "\n")
    
    # --- Region 2: Between shield and conductor (R_p < r < R) ---
    print("In the region R_p < r < R:")
    
    # H = H_r * i_r + H_theta * i_theta
    # H_r = - (2*M_0 / 3) * [ (R_p/R)^3 - (R_p/r)^3 ] * cos(theta)
    # H_theta = (M_0 / 3) * [ 2*(R_p/R)^3 + (R_p/r)^3 ] * sin(theta)
    
    h_field_outside_r_coeff = "- (2 * M_0 / 3)"
    h_field_outside_r_term = "[(R_p/R)^3 - (R_p/r)^3]"
    h_field_outside_r_vec = "cos(theta) * i_r"
    
    h_field_outside_theta_coeff = "(M_0 / 3)"
    h_field_outside_theta_term = "[2*(R_p/R)^3 + (R_p/r)^3]"
    h_field_outside_theta_vec = "sin(theta) * i_theta"

    print("The radial component H_r is:")
    print(f"  H_r = {h_field_outside_r_coeff} * {h_field_outside_r_term} * {h_field_outside_r_vec}")
    print("   Key numerical coefficients are -2/3.")
    
    print("\nThe tangential component H_theta is:")
    print(f"  H_theta = {h_field_outside_theta_coeff} * {h_field_outside_theta_term} * {h_field_outside_theta_vec}")
    print("   Key numerical coefficients are 1/3 and 2.")

# Run the display function
display_magnetic_field_solution()