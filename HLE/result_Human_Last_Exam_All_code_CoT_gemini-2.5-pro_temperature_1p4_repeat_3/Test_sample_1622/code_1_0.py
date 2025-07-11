def generate_formula():
    """
    This function generates and prints the formula for P(n).
    L represents ln(n).
    """
    
    # Numerator and denominator for the n^-2 term
    p2_num = "3*L**2 - 2*L + 2"
    p2_den = "24*n**2"
    
    # Numerator and denominator for the n^-3 term
    p3_num = "L**3 - 2*L**2 + 2*L"
    p3_den = "48*n**3"

    # Combine terms for the final formula
    formula = f"({p2_num}) / ({p2_den}) + ({p3_num}) / ({p3_den})"
    
    print("The formula for P(n) is:")
    print(f"P(n) = {formula}")

generate_formula()