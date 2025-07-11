def solve():
    """
    This function prints the derived formula for P(n).
    L represents ln(n).
    """
    p2_numerator = "3*L**2 + 2*L - 2"
    p2_denominator = "24*n**2"
    
    p3_numerator = "L**3 + 2*L**2 - 2*L"
    p3_denominator = "48*n**3"

    formula_p2 = f"({p2_numerator})/({p2_denominator})"
    formula_p3 = f"({p3_numerator})/({p3_denominator})"
    
    final_formula = f"P(n) = {formula_p2} + {formula_p3}"
    
    # For a slightly cleaner mathematical look, we'll format it without asterisks.
    pretty_formula = "(3L^2 + 2L - 2)/(24n^2) + (L^3 + 2L^2 - 2L)/(48n^3)"
    
    print(pretty_formula)

solve()