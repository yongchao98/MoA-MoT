def print_magnetic_field_solution():
    """
    This function prints the final symbolic expressions for the magnetic field H
    in the two specified regions of the substation shielding problem.
    The expressions correspond to the correct answer choice.
    """

    print("The magnetic field H(r, theta) is determined as follows:")
    print("-" * 60)

    # --- Region 1: 0 < r < R_p (within the magnetized shield) ---
    print("In the region 0 < r < R_p:")
    
    # The field is uniform and can be written in vector form.
    # The numbers in the equation are 2 and 3.
    h1_expression = "H = M_0 * ((2 * R_p**3 + R**3) / (3 * R**3)) * (-cos(theta) * i_r + sin(theta) * i_theta)"
    print(h1_expression)
    print("Numbers in this equation: 2, 3")

    print("\n" + "="*60 + "\n")

    # --- Region 2: R_p < r < R (between shield and conductor) ---
    print("In the region R_p < r < R:")

    # Radial component H_r
    # The numbers in this component's equation are -2, 3, 3, 3.
    hr2_expression = "H_r = (-2 * M_0 / 3) * ((R_p/R)**3 - (R_p/r)**3) * cos(theta)"
    print(hr2_expression)

    # Theta component H_theta
    # The numbers in this component's equation are 1, 3, 2, 3, 3.
    htheta2_expression = "H_theta = (M_0 / 3) * (2 * (R_p/R)**3 + (R_p/r)**3) * sin(theta)"
    print(htheta2_expression)
    print("Numbers in the r-component equation: -2, 3")
    print("Numbers in the theta-component equation: 1, 3, 2")


# Execute the function to display the solution
print_magnetic_field_solution()