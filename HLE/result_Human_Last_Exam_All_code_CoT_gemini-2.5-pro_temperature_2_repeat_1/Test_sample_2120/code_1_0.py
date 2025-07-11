import numpy as np

def solve_complex_average():
    """
    Calculates the average value of the complex coordinates z where the
    inverse of the gradient of S(z) approaches zero.

    This corresponds to the average of the poles of the function S(z).
    The poles of S(z) are the union of the poles of E(z) and B(z).
    """

    # Part 1: Poles of E(z)
    # Based on the analysis of the functional equation for E(z), the singularities
    # are located at z = 1 and z = 2.
    poles_E = [1, 2]
    num_poles_E = len(poles_E)
    sum_poles_E = sum(poles_E)
    print(f"Poles of E(z) are located at: {poles_E}")
    print(f"Number of poles from E(z): {num_poles_E}")
    print(f"Sum of poles from E(z): {sum_poles_E}\n")

    # Part 2: Poles of B(z)
    # The poles of B(z) are given by the roots of two polynomials, P(z) and Q(z).
    # P(z) = 4*z^4 - z^3 + z^2 + 1
    # Q(z) = z^4 + z^2 - z + 4
    
    # Coefficients of P(z) = a4*z^4 + a3*z^3 + a2*z^2 + a1*z + a0
    poly_P_coeffs = {'a4': 4, 'a3': -1, 'a2': 1, 'a1': 0, 'a0': 1}
    # Coefficients of Q(z) = b4*z^4 + b3*z^3 + b2*z^2 + b1*z + b0
    poly_Q_coeffs = {'a4': 1, 'a3': 0, 'a2': 1, 'a1': -1, 'a0': 4}

    # According to Vieta's formulas, the sum of roots of a polynomial
    # a_n*z^n + a_{n-1}*z^{n-1} + ... + a_0 is -a_{n-1}/a_n.

    # For P(z)
    num_poles_P = 4
    sum_poles_P = -poly_P_coeffs['a3'] / poly_P_coeffs['a4']
    print(f"Polynomial P(z) is {poly_P_coeffs['a4']}z^4 + {poly_P_coeffs['a3']}z^3 + {poly_P_coeffs['a2']}z^2 + {poly_P_coeffs['a1']}z + {poly_P_coeffs['a0']}")
    print(f"Sum of roots of P(z): -({poly_P_coeffs['a3']}) / {poly_P_coeffs['a4']} = {sum_poles_P}")

    # For Q(z)
    num_poles_Q = 4
    sum_poles_Q = -poly_Q_coeffs['a3'] / poly_Q_coeffs['a4']
    print(f"Polynomial Q(z) is {poly_Q_coeffs['a4']}z^4 + {poly_Q_coeffs['a3']}z^3 + {poly_Q_coeffs['a2']}z^2 + {poly_Q_coeffs['a1']}z + {poly_Q_coeffs['a0']}")
    print(f"Sum of roots of Q(z): -({poly_Q_coeffs['a3']}) / {poly_Q_coeffs['a4']} = {sum_poles_Q}\n")
    
    # Total poles from B(z)
    num_poles_B = num_poles_P + num_poles_Q
    sum_poles_B = sum_poles_P + sum_poles_Q
    print(f"Total number of poles from B(z): {num_poles_B}")
    print(f"Total sum of poles from B(z): {sum_poles_B}\n")

    # Part 3: Average of all poles
    total_num_poles = num_poles_E + num_poles_B
    total_sum_poles = sum_poles_E + sum_poles_B

    average_value = total_sum_poles / total_num_poles

    print("Final Calculation:")
    print(f"Total number of poles: {num_poles_E} + {num_poles_B} = {total_num_poles}")
    print(f"Total sum of poles: {sum_poles_E} + {sum_poles_B} = {total_sum_poles}")
    print(f"Average value = ({sum_poles_E} + {sum_poles_P} + {sum_poles_Q}) / ({num_poles_E} + {num_poles_P} + {num_poles_Q})")
    print(f"Average value = {total_sum_poles} / {total_num_poles} = {average_value}")

    return average_value

# Run the solver and print the final answer
final_average = solve_complex_average()
#<<<13/40>>>
print(f"<<<{final_average}>>>")