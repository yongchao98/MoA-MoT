def print_formula():
    """
    This function prints the derived formula for P(n).
    The formula is split into two parts for clarity.
    """
    # Using L as a placeholder for ln(n)
    P_n_term1_numerator = "3*L**2 - 2*L + 2"
    P_n_term1_denominator = "24*n**2"
    
    P_n_term2_numerator = "L**3 - 2*L**2 + 2*L"
    P_n_term2_denominator = "48*n**3"

    print("The formula for P(n) is:")
    print(f"P(n) = ({P_n_term1_numerator}) / ({P_n_term1_denominator}) + ({P_n_term2_numerator}) / ({P_n_term2_denominator})")

print_formula()