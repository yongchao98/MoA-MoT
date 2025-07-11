def solve_average_coordinate():
    """
    Calculates the average value of the complex coordinates z where the
    inverse of the gradient of S(z) approaches zero.
    """
    # From the analysis of the B(z) field equation, the singularities
    # are the roots of 4*z^4 - z^3 + z^2 + 1 = 0.
    # By Vieta's formulas, the sum of these four roots is -(-1)/4.
    sum_of_roots_B = 0.25
    num_roots_B = 4

    # From the analysis of the E(z) field equation, the singularities
    # are inferred to be at z=1 and z=2.
    root_E1 = 1
    root_E2 = 2
    num_roots_E = 2
    
    # The total sum of the coordinates is the sum of all identified singularities.
    total_sum = sum_of_roots_B + root_E1 + root_E2

    # The total number of singular coordinates.
    total_num_points = num_roots_B + num_roots_E

    # The average value is the total sum divided by the total number of points.
    average_value = total_sum / total_num_points

    print("The final calculation is based on the average of the singular points of the fields.")
    print("The equation for the average is: (sum_of_B_singularities + sum_of_E_singularities) / total_number_of_singularities")
    print(f"Sum of the 4 singular points from B(z): {sum_of_roots_B}")
    print(f"First singular point from E(z): {root_E1}")
    print(f"Second singular point from E(z): {root_E2}")
    print(f"Total number of singular points: {total_num_points}")
    print(f"The final equation is: ({sum_of_roots_B} + {root_E1} + {root_E2}) / {total_num_points}")
    print(f"Result of the calculation: {average_value}")

solve_average_coordinate()