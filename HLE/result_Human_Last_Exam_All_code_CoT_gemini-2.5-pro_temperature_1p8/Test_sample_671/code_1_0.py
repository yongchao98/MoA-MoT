def print_formula():
    """
    Prints the derived formula for the sequence of polynomials f_n(p).
    The numbers in the formula are explicitly handled as variables to meet the
    output requirements.
    """
    
    one = 1
    two = 2
    
    # The formula is f_n(p) = (p^n - (1-p)^n) / (2p - 1)
    
    formula_str = f"f_n(p) = (p^n - ({one} - p)^n) / ({two}*p - {one})"
    
    print("The simple formula for f_n(p) is:")
    print(formula_str)

print_formula()