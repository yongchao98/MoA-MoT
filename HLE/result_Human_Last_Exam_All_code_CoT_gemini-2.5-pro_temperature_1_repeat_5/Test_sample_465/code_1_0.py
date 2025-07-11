def solve():
    """
    This function determines the true statements about chromatic and orbital chromatic roots.
    A. True: The largest real orbital chromatic root is bounded by the largest real chromatic root.
    B. True: Chromatic roots can be complex numbers (e.g., for the cycle graph C4).
    C. True: Real chromatic roots can be negative, although examples are complex.
    D. True: Real chromatic roots can be non-integers (e.g., related to the golden ratio).
    E. False: Chromatic polynomials are known to have no roots in the interval (0, 1).
    """
    
    # The true statements are A, B, C, and D.
    # The required output is a sorted string of the letters corresponding to the true statements.
    true_statements = ['A', 'B', 'C', 'D']
    
    # Sort the letters alphabetically
    true_statements.sort()
    
    # Join them into a single string
    result = "".join(true_statements)
    
    print(result)

solve()