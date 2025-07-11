def solve_chromatic_roots_properties():
    """
    This function evaluates the truth of given statements about chromatic
    and orbital chromatic roots and constructs the final answer string.
    """
    # Step 1: Analyze each statement based on known theorems in graph theory.
    # Statement A: True. Known as the Interlacing Theorem for Orbital Chromatic Polynomials.
    # Statement B: True. Many graphs have chromatic polynomials with non-real complex roots.
    # Statement C: False. It is known that P(G, k) != 0 for k < 0.
    # Statement D: True. For example, φ^2 = (3+sqrt(5))/2 is a chromatic root for some graphs.
    # Statement E: False. It is a theorem that there are no chromatic roots in the interval (0, 1).

    statements_truth_values = {
        'A': True,
        'B': True,
        'C': False,
        'D': True,
        'E': False
    }

    # Step 2: Collect the letters of the true statements.
    true_statements = []
    for statement, is_true in statements_truth_values.items():
        if is_true:
            true_statements.append(statement)

    # Step 3: Sort the letters alphabetically.
    true_statements.sort()

    # Step 4: Join the letters into a single string.
    result = "".join(true_statements)

    # Step 5: If no statements were true, the answer is "0".
    if not result:
        result = "0"
        
    print(f"The true statements are: {', '.join(true_statements)}")
    print(f"The final answer as a sorted string is: {result}")

solve_chromatic_roots_properties()