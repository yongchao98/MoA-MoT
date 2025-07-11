def solve():
    """
    This function determines which of the given statements about chromatic
    and orbital chromatic roots are true, and prints the sorted result.
    
    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (False)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1. (False)
    
    The true statements are A, B, and D.
    """
    
    true_statements = ['A', 'B', 'D']
    
    # Sort the list of true statements alphabetically
    true_statements.sort()
    
    # Join the characters to form the final string
    answer = "".join(true_statements)
    
    print(answer)

solve()