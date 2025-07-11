def solve_quiz():
    """
    Analyzes the given statements about chromatic and orbital chromatic roots
    and prints the sorted string of true statements.

    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (True)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1. (False)
    """
    
    true_statements = ['A', 'B', 'C', 'D']
    
    # Sort the letters alphabetically
    true_statements.sort()
    
    # Join them into a single string
    result = "".join(true_statements)
    
    print(result)

solve_quiz()