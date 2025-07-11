def print_p_formula():
    """
    This function prints the formula for P(n) with its numerical coefficients.
    L is defined as ln(n).
    """

    # Coefficients for the first term of P(n), with n^2 in the denominator.
    # The numerator is c22*L^2 + c21*L + c20
    c22 = 3
    c21 = -2
    c20 = 2
    d2 = 24

    # Coefficients for the second term of P(n), with n^3 in the denominator.
    # The numerator is c33*L^3 + c32*L^2 + c31*L
    c33 = 1
    c32 = -2
    c31 = 2
    d3 = 48
    
    # Construct and print the formula string
    formula = f"({c22}*L^2 - {abs(c21)}*L + {c20})/({d2}*n^2) + ({c33}*L^3 - {abs(c32)}*L^2 + {c31}*L)/({d3}*n^3)"
    print(formula)

print_p_formula()