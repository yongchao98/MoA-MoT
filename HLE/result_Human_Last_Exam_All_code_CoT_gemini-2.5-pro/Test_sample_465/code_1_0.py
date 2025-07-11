def solve_chromatic_roots_question():
    """
    This function determines the true statements about chromatic and orbital chromatic roots
    and prints the answer in the required format.

    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (True)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1. (False)

    The true statements are A, B, C, and D.
    The sorted string of these letters is "ABCD".
    """
    
    # The letters corresponding to the true statements.
    true_statements = ['A', 'B', 'C', 'D']
    
    # Sort the letters alphabetically and join them into a single string.
    answer = "".join(sorted(true_statements))
    
    # Print the final answer.
    print(answer)

solve_chromatic_roots_question()