def solve_chromatic_roots_properties():
    """
    This function analyzes the properties of chromatic and orbital chromatic roots.

    The statements are:
    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (True)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1 (i.e., r in (0, 1)). (False)

    The true statements are A, B, C, D.
    The required output is a sorted string of the letters corresponding to the true statements.
    """
    
    true_statements = ['A', 'B', 'C', 'D']
    
    # Sort the letters alphabetically
    sorted_answer = "".join(sorted(true_statements))
    
    print(sorted_answer)

solve_chromatic_roots_properties()