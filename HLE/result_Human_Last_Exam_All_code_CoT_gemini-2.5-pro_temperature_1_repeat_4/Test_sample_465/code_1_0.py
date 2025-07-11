def solve_chromatic_roots_properties():
    """
    Analyzes several statements about chromatic and orbital chromatic roots
    and determines which are true.
    """

    # A dictionary to hold the truth value of each statement.
    # The value is a tuple: (is_true, reason).
    statements_analysis = {
        'A': (False, "False. Counterexample: For the null graph N_n, the largest chromatic root is 0, while the largest orbital chromatic root (under the symmetric group) is n-1."),
        'B': (True, "True. The chromatic polynomial for the cycle graph C4 has non-real roots 3/2 +/- i*sqrt(3)/2."),
        'C': (True, "True. It has been shown that certain non-planar graphs have negative chromatic roots (e.g., a root near -2.43 was found by Royle and Sokal)."),
        'D': (True, "True. For a family of planar graphs, the golden ratio squared (phi^2 approx 2.618) is a chromatic root."),
        'E': (False, "False. It is a proven theorem (Jackson, 1993) that no chromatic polynomial has a root in the interval (0, 1).")
    }

    true_statements = []
    for statement, (is_true, reason) in sorted(statements_analysis.items()):
        if is_true:
            true_statements.append(statement)

    # The final answer is a sorted string of the letters of true statements.
    answer = "".join(true_statements)

    if not answer:
        answer = "0"

    print(answer)

solve_chromatic_roots_properties()