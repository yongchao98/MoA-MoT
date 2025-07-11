def solve_chromatic_roots_problem():
    """
    This function determines the true statements about chromatic and orbital chromatic roots
    from the given list and prints the answer in the specified format.

    The statements are:
    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (False - this is a famous open problem)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1 (i.e., r in (0, 1)). (False)

    The true statements are A, B, and D.
    """
    
    # List of true statements
    true_statements = ['A', 'B', 'D']
    
    # The problem asks for the answer as a sorted string without any delimiters.
    true_statements.sort()
    
    # Join the sorted list into a single string.
    answer = "".join(true_statements)
    
    # Print the final answer.
    print(answer)

solve_chromatic_roots_problem()