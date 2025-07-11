def display_magnetic_field_equations():
    """
    Displays the equations for the magnetic field H inside and outside the shield.
    The equations correspond to the correct answer choice.
    The derivation involves solving Laplace's equation for the magnetic scalar potential
    with the appropriate boundary conditions for a magnetized sphere inside a
    perfectly conducting cavity.
    """
    # Define the equations as strings for clear printing
    region_1_condition = "In the region 0 < r < R_p:"
    eq_1 = "H = M_0 * ((2*R_p**3 + R**3) / (3*R**3)) * (-cos(theta) * i_r + sin(theta) * i_theta)"

    region_2_condition = "In the region R_p < r < R:"
    # Breaking down the long equation for readability
    eq_2_part1 = "H = - (2*M_0 / 3) * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta) * i_r"
    eq_2_part2 = "    + (M_0 / 3) * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta) * i_theta"

    # Print the final result
    print("The magnetic field H(r, theta) is determined in two regions:")
    print("----------------------------------------------------------------------")
    print(region_1_condition)
    print("  " + eq_1)
    print("----------------------------------------------------------------------")
    print(region_2_condition)
    print("  " + eq_2_part1)
    print("  " + eq_2_part2)
    print("----------------------------------------------------------------------")
    print("\nThis result matches the expressions given in Answer Choice B.")

# Execute the function to display the answer
display_magnetic_field_equations()