def solve_chromatic_roots_properties():
    """
    This function determines which statements about chromatic and orbital chromatic roots are true
    and prints the sorted combination of the corresponding letters.

    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (True)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (False)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1. (False)
    """

    # The letters corresponding to the true statements.
    true_statements = ['A', 'B', 'D']

    # Sort the letters alphabetically.
    true_statements.sort()

    # Join them into a single string.
    result = "".join(true_statements)

    print(result)

solve_chromatic_roots_properties()