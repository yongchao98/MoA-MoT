def solve_chromatic_roots_question():
    """
    Analyzes statements about chromatic and orbital chromatic roots and determines the true ones.
    """
    
    # A dictionary to hold the truth value for each statement based on established theorems in graph theory.
    # A. True: A theorem by B. Jackson (2008) confirms this.
    # B. True: A result by A. Sokal (2001) shows chromatic roots are dense in the complex plane.
    # C. False: A theorem by B. Jackson (1993) shows no chromatic roots are in (-inf, 0).
    # D. True: W. T. Tutte showed that values like the golden ratio plus one can be chromatic roots.
    # E. False: A theorem by B. Jackson (1993) shows no chromatic roots are in (0, 1).
    statements = {
        'A': True,
        'B': True,
        'C': False,
        'D': True,
        'E': False
    }
    
    # Collect the letters of the true statements.
    true_statements = []
    for letter, is_true in statements.items():
        if is_true:
            true_statements.append(letter)
            
    # Sort the letters alphabetically.
    true_statements.sort()
    
    # Join the letters into a single string.
    answer = "".join(true_statements)
    
    print(answer)

solve_chromatic_roots_question()