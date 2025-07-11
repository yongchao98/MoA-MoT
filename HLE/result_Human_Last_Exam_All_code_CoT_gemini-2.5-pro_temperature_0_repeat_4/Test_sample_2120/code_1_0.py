def solve_average_z():
    """
    Calculates the average value of the complex coordinates z where the
    inverse of the gradient of S(z) approaches zero.
    """
    # Step 1: Identify the poles of E(z).
    # From the functional equation for E(z), the poles are located at z=1 and z=2.
    pole_E1 = 1
    pole_E2 = 2
    
    # Step 2: Identify the poles of B(z).
    # The poles of B(z) are the roots of the denominator of the right-hand side
    # of its governing equation: 4*z^4 - z^3 + z^2 + 1 = 0.
    # The coefficients are a=4, b=-1, c=1, d=0, e=1.
    a = 4
    b = -1
    
    # Step 3: Calculate the sum of the poles of B(z) using Vieta's formulas.
    # The sum of the roots of a quartic polynomial ax^4+bx^3+...=0 is -b/a.
    sum_poles_B = -b / a
    
    # Step 4: Calculate the sum of all poles.
    sum_all_poles = pole_E1 + pole_E2 + sum_poles_B
    
    # Step 5: Calculate the total number of poles.
    num_poles_E = 2
    num_poles_B = 4  # A quartic polynomial has 4 roots.
    total_num_poles = num_poles_E + num_poles_B
    
    # Step 6: Compute the average value.
    average_z = sum_all_poles / total_num_poles
    
    # Print the calculation steps and the final result.
    print("The average value of the complex coordinates z is calculated by averaging the poles of E(z) and B(z).")
    print(f"Poles from E(z): {pole_E1}, {pole_E2}")
    print(f"The sum of poles from B(z) (roots of 4z^4 - z^3 + z^2 + 1 = 0) is: {sum_poles_B}")
    print(f"The average is calculated as: ({pole_E1} + {pole_E2} + {sum_poles_B}) / {total_num_poles}")
    print(f"The final average value is: {average_z}")

solve_average_z()