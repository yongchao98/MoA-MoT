def solve_chromatic_roots_properties():
    """
    Analyzes statements about chromatic and orbital chromatic roots and determines which are true.
    """
    # Statement A: Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    # Statement B: Chromatic roots may not be real. (True)
    # Statement C: Real chromatic roots may take on negative values. (False)
    # Statement D: Real chromatic roots may take on non-integer values. (True)
    # Statement E: Chromatic polynomials may have roots between 0 and 1. (False)

    true_statements = ['A', 'B', 'D']
    
    # Sort the letters alphabetically
    true_statements.sort()
    
    # Join them into a single string
    answer = "".join(true_statements)
    
    print(answer)

solve_chromatic_roots_properties()