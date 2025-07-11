import sys

def solve_complex_average():
    """
    Calculates the average value of complex coordinates z where the inverse
    of the gradient of S(z) approaches zero.
    """
    # Step 1: Problem interpretation
    # The coordinates 'z' where the inverse of the gradient of S(z) approaches zero
    # are the poles of the derivative S'(z), as |S'(z)| -> infinity.
    # The poles of S'(z) are located at the poles of the fields E(z) and B(z).

    # Step 2: Analyze the equation for B(z) to find its poles.
    # The equation for B(z) has a rational function on the right-hand side:
    # RHS = (z + z^2) / (4z^4 - z^3 + z^2 + 1)
    # The poles of B(z) are assumed to be the roots of the denominator polynomial.
    # P(z) = 4z^4 - z^3 + z^2 + 1 = 0
    poly_coeffs = {'a4': 4, 'a3': -1, 'a2': 1, 'a1': 0, 'a0': 1}
    num_b_poles = 4

    # Using Vieta's formulas, the sum of the roots of a polynomial
    # a_n*z^n + ... + a_0 = 0 is -a_{n-1}/a_n.
    sum_b_poles = -poly_coeffs['a3'] / poly_coeffs['a4']

    # Step 3: Analyze the equation for E(z) to find its poles.
    # The functional equation for E(z) involves arguments f1(z) = (z-3)/(z-2)
    # and f2(z) = (2z-3)/(z-1).
    # These arguments have poles at z=2 and z=1, respectively.
    # This structure implies that E(z) must also have poles at these locations
    # to satisfy the equation.
    pole_e1 = 1
    pole_e2 = 2
    num_e_poles = 2

    # Step 4: Calculate the average of all identified pole coordinates.
    total_num_poles = num_b_poles + num_e_poles
    sum_of_all_poles = sum_b_poles + pole_e1 + pole_e2
    average_z = sum_of_all_poles / total_num_poles

    # Step 5: Output the detailed reasoning and the final calculation.
    print("This script calculates the average of complex coordinates 'z' where the gradient of S(z) has poles.")
    print("This is equivalent to finding the average of the poles of the E(z) and B(z) fields.")
    
    print("\n--- Poles from B(z) ---")
    print(f"The poles of B(z) are the roots of the polynomial: {poly_coeffs['a4']}z^4 + ({poly_coeffs['a3']})z^3 + {poly_coeffs['a2']}z^2 + {poly_coeffs['a0']} = 0.")
    print(f"Number of poles from B(z): {num_b_poles}")
    print("The sum of these roots is calculated using Vieta's formulas as -a3/a4.")
    print(f"Sum of poles from B(z) = -({poly_coeffs['a3']}) / {poly_coeffs['a4']} = {sum_b_poles}")

    print("\n--- Poles from E(z) ---")
    print("The functional equation for E(z) implies singularities based on its arguments' denominators.")
    print(f"Poles for E(z) are identified at z = {pole_e1} and z = {pole_e2}.")
    print(f"Number of poles from E(z): {num_e_poles}")

    print("\n--- Average Calculation ---")
    print(f"Total number of poles = {num_b_poles} + {num_e_poles} = {total_num_poles}")
    print(f"Total sum of all poles = {sum_b_poles} + {pole_e1} + {pole_e2} = {sum_of_all_poles}")
    print("The average value of z is the total sum of poles divided by the total number of poles.")
    
    # Final equation as requested
    print("\nFinal Equation:")
    print(f"({sum_b_poles} + {pole_e1} + {pole_e2}) / ({num_b_poles} + {num_e_poles}) = {average_z}")

solve_complex_average()