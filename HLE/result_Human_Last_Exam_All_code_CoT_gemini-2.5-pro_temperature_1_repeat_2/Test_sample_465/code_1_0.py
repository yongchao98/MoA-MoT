def solve_chromatic_roots_quiz():
    """
    Determines which statements about chromatic roots are true based on known mathematical theorems.
    """

    # A dictionary to hold the analysis for each statement.
    # The truth value is determined by established results in graph theory.
    statements_analysis = {
        'A': {
            'is_true': True,
            'reason': "Theorem by Cameron, Morgan, and Jackson shows real orbital chromatic roots are in [0, ρ(G)]."
        },
        'B': {
            'is_true': True,
            'reason': "Chromatic roots can be complex. The set of roots is dense in the complex plane (Sokal's theorem)."
        },
        'C': {
            'is_true': False,
            'reason': "Real chromatic roots are non-negative (Jackson's theorem)."
        },
        'D': {
            'is_true': True,
            'reason': "Non-integer roots exist, e.g., Beraha numbers like (3+sqrt(5))/2."
        },
        'E': {
            'is_true': False,
            'reason': "Chromatic polynomials have no roots in the interval (0, 1) (Jackson's theorem)."
        }
    }

    # Collect the identifiers of the true statements.
    true_statements = []
    for statement_id, analysis in statements_analysis.items():
        if analysis['is_true']:
            true_statements.append(statement_id)

    # Sort the identifiers alphabetically and join them into a string.
    # If no statements are true, the result should be "0".
    if not true_statements:
        final_answer = "0"
    else:
        # Sort for a consistent, ordered output.
        true_statements.sort()
        final_answer = "".join(true_statements)

    print(final_answer)

solve_chromatic_roots_quiz()