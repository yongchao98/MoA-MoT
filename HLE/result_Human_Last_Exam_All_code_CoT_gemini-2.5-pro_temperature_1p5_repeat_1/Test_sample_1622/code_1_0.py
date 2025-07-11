def generate_formula():
    """
    This function generates the formula for P(n) based on the derived coefficients.
    """
    # Coefficients for the first term of P(n), which is of order 1/n^2.
    # The numerator is 3*L**2 + 2*L - 2.
    p2_coeff_L2 = 3
    p2_coeff_L1 = 2
    p2_coeff_L0 = -2
    p2_denominator = 24

    # Coefficients for the second term of P(n), which is of order 1/n^3.
    # The numerator is L**3 + 2*L**2 - 2*L.
    p3_coeff_L3 = 1
    p3_coeff_L2 = 2
    p3_coeff_L1 = -2
    p3_denominator = 48

    # Construct the string for the first term
    term1_numerator = f"({p2_coeff_L2}*L**2 + {p2_coeff_L1}*L - {abs(p2_coeff_L0)})"
    term1_denominator = f"({p2_denominator}*n**2)"
    term1 = f"{term1_numerator}/{term1_denominator}"
    
    # Construct the string for the second term
    # Handling the case where the leading coefficient is 1.
    term2_numerator = f"L**3 + {p3_coeff_L2}*L**2 - {abs(p3_coeff_L1)}*L"
    term2_numerator_parenthesized = f"({term2_numerator})"
    term2_denominator = f"({p3_denominator}*n**3)"
    term2 = f"{term2_numerator_parenthesized}/{term2_denominator}"
    
    # Combine the terms to form the full formula for P(n)
    final_formula = f"{term1} + {term2}"

    print(final_formula)

generate_formula()