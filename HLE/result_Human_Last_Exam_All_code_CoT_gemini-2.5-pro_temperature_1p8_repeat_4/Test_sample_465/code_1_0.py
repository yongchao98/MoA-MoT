def solve_chromatic_roots_question():
    """
    This function determines which statements about chromatic and orbital chromatic roots are true.

    The statements are:
    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (False)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (True)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1 (i.e., r in (0, 1)). (False)

    The true statements are B, C, and D.
    The function returns these letters as a sorted string.
    """
    true_statements = ['B', 'C', 'D']
    
    # Sort the list of characters alphabetically
    true_statements.sort()
    
    # Join the characters to form the final answer string
    answer = "".join(true_statements)
    
    print(answer)

solve_chromatic_roots_question()