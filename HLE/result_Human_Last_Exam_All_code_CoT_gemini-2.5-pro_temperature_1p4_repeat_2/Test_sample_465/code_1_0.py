def solve_chromatic_roots_properties():
    """
    Analyzes statements about chromatic and orbital chromatic roots and prints the conclusion.

    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (False)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1. (False)

    The true statements are A, B, and D.
    """
    true_statements = ['A', 'B', 'D']
    
    # Sort the letters alphabetically and join them into a string.
    answer = "".join(sorted(true_statements))
    
    print(answer)

solve_chromatic_roots_properties()