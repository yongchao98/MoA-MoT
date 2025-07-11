import numpy as np

def solve_coordinates_average():
    """
    Calculates the average value of the complex coordinates z based on the
    singularities derived from the provided B-field equation.

    The plan is to find the average of the roots of the denominator polynomial
    of the source term in the B-field equation, as these roots are assumed
    to be the singularities of the system.
    """

    # The denominator polynomial from the B-field equation is P(z) = 4z^4 - z^3 + z^2 + 1.
    # The coefficients are [4, -1, 1, 0, 1] for powers z^4, z^3, z^2, z^1, z^0.
    coeffs = [4, -1, 1, 0, 1]
    
    # According to Vieta's formulas, for a polynomial of degree n,
    # P(z) = a_n*z^n + a_{n-1}*z^{n-1} + ... + a_0
    # the sum of the roots is -a_{n-1} / a_n.
    
    n = len(coeffs) - 1
    a_n = coeffs[0]
    a_n_minus_1 = coeffs[1]
    
    # Calculate the sum of the roots.
    sum_of_roots = -a_n_minus_1 / a_n
    
    # The number of roots is the degree of the polynomial.
    number_of_roots = n
    
    # Calculate the average of the roots.
    average_of_roots = sum_of_roots / number_of_roots
    
    # The final equation is Average = (Sum of Roots) / (Number of Roots)
    # Let's print the numbers in the final equation.
    
    print(f"The denominator polynomial is P(z) = {a_n}z^4 + ({a_n_minus_1})z^3 + ...")
    print(f"The sum of the roots is calculated as -({a_n_minus_1}) / {a_n} = {sum_of_roots}")
    print(f"The number of roots is the degree of the polynomial, which is {number_of_roots}")
    
    print("\nThe final equation for the average is:")
    # We use round to avoid floating point representation issues for 0.25
    print(f"Average = {round(sum_of_roots, 4)} / {number_of_roots}")
    
    print(f"\nThe average value of the complex coordinates is: {average_of_roots}")

solve_coordinates_average()