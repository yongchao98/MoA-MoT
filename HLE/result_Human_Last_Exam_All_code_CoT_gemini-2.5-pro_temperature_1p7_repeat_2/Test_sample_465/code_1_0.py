def solve_chromatic_roots_properties():
    """
    Analyzes the properties of chromatic and orbital chromatic roots and prints the final answer.

    The statements are:
    A. Real orbital chromatic roots are bounded by the greatest real chromatic root. (False)
    B. Chromatic roots may not be real. (True)
    C. Real chromatic roots may take on negative values. (True)
    D. Real chromatic roots may take on non-integer values. (True)
    E. Chromatic polynomials may have roots between 0 and 1. (False)
    """

    # Dictionary to store the truth value of each statement
    statements_truth = {
        'A': False,
        'B': True,
        'C': True,
        'D': True,
        'E': False
    }

    true_statements = []
    for statement, is_true in statements_truth.items():
        if is_true:
            true_statements.append(statement)

    # Sort the letters alphabetically
    true_statements.sort()

    # Join to form the final string, or "0" if none are true
    if not true_statements:
        final_answer = "0"
    else:
        final_answer = "".join(true_statements)

    print(final_answer)

solve_chromatic_roots_properties()