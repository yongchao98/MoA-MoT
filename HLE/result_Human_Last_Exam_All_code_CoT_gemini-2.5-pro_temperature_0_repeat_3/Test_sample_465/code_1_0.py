def solve_chromatic_roots_problem():
    """
    Analyzes statements about chromatic and orbital chromatic roots
    and determines which are true.
    """

    statements = {
        'A': {
            'text': "Real orbital chromatic roots are bounded by the greatest real chromatic root.",
            'is_true': True,
            'justification': "This is a known result. The largest real root of the orbital chromatic polynomial P(G, Γ) is less than or equal to the largest real root of the chromatic polynomial P(G)."
        },
        'B': {
            'text': "Chromatic roots may not be real.",
            'is_true': True,
            'justification': "This is true. Chromatic polynomials can have complex roots. The Petersen graph is a well-known example."
        },
        'C': {
            'text': "Real chromatic roots may take on negative values.",
            'is_true': True,
            'justification': "This is true. It has been shown that the set of real chromatic roots is dense in the interval (-∞, 32/27], which includes negative values."
        },
        'D': {
            'text': "Real chromatic roots may take on non-integer values.",
            'is_true': True,
            'justification': "This is true. W. T. Tutte showed that for a certain family of planar graphs, the square of the golden ratio, (3 + sqrt(5))/2, is a chromatic root. This value is not an integer."
        },
        'E': {
            'text': "Chromatic polynomials may have roots between 0 and 1 (i.e., r ∈ (0, 1)).",
            'is_true': False,
            'justification': "This is false. A theorem by Jackson (1993) proves that chromatic polynomials have no roots in the interval (0, 1)."
        }
    }

    true_statements = []
    for statement_letter, details in statements.items():
        if details['is_true']:
            true_statements.append(statement_letter)

    # Sort the letters alphabetically
    true_statements.sort()

    # Join them into a single string
    final_answer = "".join(true_statements)

    print(final_answer)

solve_chromatic_roots_problem()