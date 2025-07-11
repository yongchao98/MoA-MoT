def find_true_statements():
    """
    This function determines which statements about chromatic and orbital chromatic roots are true
    based on established theorems in graph theory.
    """
    # A. Real orbital chromatic roots are bounded by the greatest real chromatic root.
    # Status: True. This is a theorem by Cameron and Cohen (1992).
    statement_A = True

    # B. Chromatic roots may not be real.
    # Status: True. Chromatic polynomials can have complex roots. For example, the diamond graph.
    statement_B = True

    # C. Real chromatic roots may take on negative values.
    # Status: False. The alternating signs of the chromatic polynomial's coefficients
    # ensure it has no roots for k < 0.
    statement_C = False

    # D. Real chromatic roots may take on non-integer values.
    # Status: True. The set of real chromatic roots is dense in [32/27, infinity),
    # which contains non-integers.
    statement_D = True

    # E. Chromatic polynomials may have roots between 0 and 1.
    # Status: False. It is a known theorem that there are no chromatic roots in the interval (0, 1).
    statement_E = False

    true_statements = []
    if statement_A:
        true_statements.append("A")
    if statement_B:
        true_statements.append("B")
    if statement_C:
        true_statements.append("C")
    if statement_D:
        true_statements.append("D")
    if statement_E:
        true_statements.append("E")

    # The final answer is a sorted string of the letters of the true statements.
    true_statements.sort()
    answer = "".join(true_statements)
    print(answer)

find_true_statements()