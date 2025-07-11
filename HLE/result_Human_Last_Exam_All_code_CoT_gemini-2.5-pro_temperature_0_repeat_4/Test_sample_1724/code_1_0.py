def find_frequency_correction_term():
    """
    This function calculates the third term of the polynomial found in the
    second-order nonlinear frequency correction for the Rayleigh-Plesset equation.

    The second-order frequency correction term, omega_2, is proportional to a
    polynomial in the polytropic index, gamma. The polynomial is:
    P(gamma) = 6*gamma^2 + 5*gamma + 14
    
    This code identifies the coefficients of this polynomial and prints the
    third term.
    """

    # The polynomial in gamma is P(gamma) = c1*gamma^2 + c2*gamma + c3
    c1 = 6
    c2 = 5
    c3 = 14

    print("The polynomial part of the frequency correction is:")
    print(f"({c1})*gamma^2 + ({c2})*gamma + ({c3})")
    
    print("\nThe first term of the polynomial is 6*gamma^2.")
    print("The second term of the polynomial is 5*gamma.")
    print("The third term of the polynomial is 14.")

    # The question asks for the 3rd term of the nonlinear correction.
    # Based on the interpretation that this refers to the third term of the
    # polynomial P(gamma), the answer is the constant term.
    answer = c3
    
    print(f"\nThe numerical value of the third term is:")
    print(answer)

find_frequency_correction_term()