def print_formula():
    """
    This function prints the derived formula for P(n).
    The derivation is based on matching the asymptotic series of the function
    and its approximation.
    """
    
    # Numerator and denominator for the n^-2 term
    num_n2_c3 = 3
    num_n2_c2 = -2
    num_n2_c1 = 2
    den_n2 = 24
    
    # Numerator and denominator for the n^-3 term
    num_n3_c4 = 1
    num_n3_c3 = -2
    num_n3_c2 = 2
    den_n3 = 48
    
    # Using an f-string to construct and print the formula clearly.
    # The variable L represents ln(n).
    formula_str = (
        f"P(n) = ({num_n2_c3}*L^2 + ({num_n2_c2})*L + {num_n2_c1}) / ({den_n2}*n^2) + "
        f"({num_n3_c4}*L^3 + ({num_n3_c3})*L^2 + {num_n3_c2}*L) / ({den_n3}*n^3)"
    )
    
    print(formula_str)

print_formula()