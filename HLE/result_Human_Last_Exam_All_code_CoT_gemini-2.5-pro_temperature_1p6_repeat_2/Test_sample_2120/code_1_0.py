def solve_complex_average():
    """
    Calculates the average value of the complex coordinates z where the inverse
    of the gradient of S(z) approaches zero.
    """
    
    # According to Vieta's formulas, the sum of the roots of a polynomial
    # a_n*z^n + ... + a_0 = 0 is -a_{n-1}/a_n.
    # For the polynomial 4*z^4 - z^3 + z^2 + 1 = 0 from the B-field equation:
    # a_4 = 4, a_3 = -1.
    sum_b_poles = -(-1) / 4.0
    num_b_poles = 4
    
    # The poles from the E-field equation are at z=1 and z=2.
    pole_e1 = 1
    pole_e2 = 2
    sum_e_poles = pole_e1 + pole_e2
    num_e_poles = 2
    
    # The total sum and total number of poles
    total_sum_of_poles = sum_b_poles + sum_e_poles
    total_number_of_poles = num_b_poles + num_e_poles
    
    # Calculate the average
    average_value = total_sum_of_poles / total_number_of_poles
    
    # Output the result showing the components of the calculation.
    # The problem requests to output each number in the final equation.
    print(f"The average is calculated as the sum of all pole coordinates divided by the number of poles.")
    print(f"Sum of poles from B-field equation (roots of 4z^4-z^3+z^2+1=0): {sum_b_poles}")
    print(f"Poles from E-field equation: {pole_e1}, {pole_e2}")
    print(f"Number of poles from B-field: {num_b_poles}")
    print(f"Number of poles from E-field: {num_e_poles}")
    print("\nFinal calculation:")
    print(f"({sum_b_poles} + {pole_e1} + {pole_e2}) / ({num_b_poles} + {num_e_poles}) = {total_sum_of_poles} / {total_number_of_poles} = {average_value}")

solve_complex_average()
print(f"\n<<<{13/24}>>>")