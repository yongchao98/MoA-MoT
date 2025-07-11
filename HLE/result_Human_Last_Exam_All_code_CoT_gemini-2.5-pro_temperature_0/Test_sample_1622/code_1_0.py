def solve():
    """
    This function prints the formula for P(n).
    L is used as a shorthand for ln(n).
    """
    
    # The formula for P(n) is derived by matching the asymptotic expansions.
    # It consists of terms of order n^-2 and n^-3.
    term1_numerator = "3*L**2 - 2*L + 2"
    term1_denominator = "24*n**2"
    
    term2_numerator = "L**3 - 2*L**2 + 2*L"
    term2_denominator = "48*n**3"
    
    formula = f"P(n) = ({term1_numerator}) / ({term1_denominator}) + ({term2_numerator}) / ({term2_denominator})"
    
    print(formula)

solve()