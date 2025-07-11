def find_probability_formula():
    """
    This function prints the derived formula for P_m.
    The formula depends on the parity of m.
    """
    
    print("The probability P_m is given by the formula:")
    
    # Case for odd m
    numerator_odd = 3
    denominator_expr = "(2m+1)(4m+1)"
    print(f"If m is odd: P_m = {numerator_odd} / {denominator_expr}")
    
    # Case for even m
    numerator_even = 4
    print(f"If m is even: P_m = {numerator_even} / {denominator_expr}")
    
    print("\nBreaking down the components of the formula:")
    print(f"Numerator: {numerator_odd} (if m is odd), {numerator_even} (if m is even)")
    print("Denominator terms: (2m + 1) and (4m + 1)")

find_probability_formula()