def solve_chromatic_roots_quiz():
    """
    This function determines which statements about chromatic and orbital chromatic roots are true
    and prints the sorted result.

    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (False)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1. (False)
    """
    
    # List of true statements
    true_statements = ['A', 'B', 'D']
    
    # Sort the list alphabetically and join into a single string
    answer = "".join(sorted(true_statements))
    
    print(answer)

solve_chromatic_roots_quiz()