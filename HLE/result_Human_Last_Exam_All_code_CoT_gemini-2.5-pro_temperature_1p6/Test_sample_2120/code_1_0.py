import numpy as np

def solve_field_average():
    """
    Calculates the average value of the complex coordinates z where the
    inverse of the gradient of S(z) approaches zero.

    This is equivalent to finding the average of the poles of S'(z).
    We assume the poles are dominated by the B(z) field, as the E(z)
    equation's structure suggests it might be an entire function.

    The poles of B(z) come from the denominators of the functions h(z) and h(1/z).
    h(z) = (z+z^2) / (1+z^2-z^3+4z^4)
    Poles of h(z) are roots of P1(z) = 4z^4 - z^3 + z^2 + 1 = 0.
    Poles of h(1/z) are roots of P2(z) = z^4 + z^2 - z + 4 = 0.
    The total set of poles for B(z) are the roots of P(z) = P1(z) * P2(z).
    """

    # Coefficients of P1(z) = 4z^4 - z^3 + z^2 + 1
    # from highest degree to lowest
    p1_coeffs = [4, -1, 1, 0, 1]
    
    # Coefficients of P2(z) = z^4 + z^2 - z + 4
    # from highest degree to lowest
    p2_coeffs = [1, 0, 1, -1, 4]

    # Represent as numpy polynomial objects
    P1 = np.polynomial.Polynomial(p1_coeffs[::-1])
    P2 = np.polynomial.Polynomial(p2_coeffs[::-1])
    
    # Multiply the polynomials to get the polynomial for all poles
    P_total = P1 * P2
    
    # Get the coefficients of the resulting polynomial
    total_coeffs = P_total.coef[::-1]
    
    # The degree of the polynomial is the number of poles
    number_of_poles = P_total.degree()
    
    # According to Vieta's formulas, the sum of the roots of a polynomial
    # a_n*z^n + a_{n-1}*z^{n-1} + ... + a_0 = 0
    # is -a_{n-1} / a_n.
    # a_n is the coefficient of the highest degree term (total_coeffs[0])
    # a_{n-1} is the coefficient of the second highest degree term (total_coeffs[1])
    
    coeff_n = total_coeffs[0]
    coeff_n_minus_1 = total_coeffs[1]
    
    sum_of_poles = -coeff_n_minus_1 / coeff_n
    
    # Calculate the average value
    average_value = sum_of_poles / number_of_poles
    
    print(f"The polynomial whose roots are the poles is given by the product of P1(z) and P2(z).")
    print(f"P1(z) = {p1_coeffs[0]}z^4 + {p1_coeffs[1]}z^3 + {p1_coeffs[2]}z^2 + {p1_coeffs[3]}z + {p1_coeffs[4]}")
    print(f"P2(z) = {p2_coeffs[0]}z^4 + {p2_coeffs[1]}z^3 + {p2_coeffs[2]}z^2 + {p2_coeffs[3]}z + {p2_coeffs[4]}")
    
    print("\nCalculating the sum and number of poles:")
    print(f"Sum of poles = -({coeff_n_minus_1}) / {coeff_n} = {sum_of_poles}")
    print(f"Number of poles = {number_of_poles}")

    print("\nThe final equation for the average is:")
    print(f"Average value = Sum of poles / Number of poles")
    print(f"Average value = {sum_of_poles} / {number_of_poles} = {average_value}")

solve_field_average()