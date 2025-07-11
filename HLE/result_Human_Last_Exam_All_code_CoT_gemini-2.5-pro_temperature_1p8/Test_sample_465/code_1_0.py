def solve_chromatic_roots_properties():
    """
    Analyzes statements about chromatic and orbital chromatic roots
    and determines which are true.
    """
    
    # A dictionary to hold the truth value and a brief explanation for each statement.
    statements = {
        'A': {
            "text": "Real orbital chromatic roots are bounded by the greatest real chromatic root.",
            "is_true": True,
            "reason": "This is a known theorem by Cameron and Jackson."
        },
        'B': {
            "text": "Chromatic roots may not be real.",
            "is_true": True,
            "reason": "Chromatic roots can be complex and are known to be dense in the complex plane."
        },
        'C': {
            "text": "Real chromatic roots may take on negative values.",
            "is_true": True,
            "reason": "Read's conjecture (that they are non-negative) was disproven with counterexamples."
        },
        'D': {
            "text": "Real chromatic roots may take on non-integer values.",
            "is_true": True,
            "reason": "The Tutte-Beraha numbers, which are not all integers, are limit points of real chromatic roots."
        },
        'E': {
            "text": "Chromatic polynomials may have roots between 0 and 1 (i.e., r in (0, 1)).",
            "is_true": False,
            "reason": "A fundamental theorem by Jackson and Thomassen proves that there are no chromatic roots in the interval (0, 1)."
        }
    }
    
    # Collect the labels of the true statements.
    true_statements = []
    for label, info in statements.items():
        if info["is_true"]:
            true_statements.append(label)
            
    # Sort the labels alphabetically.
    true_statements.sort()
    
    # Join them into a single string for the final answer.
    # If no statements were true, the result would be an empty list,
    # and we would output '0'. In this case, we have true statements.
    if not true_statements:
        final_answer = "0"
    else:
        final_answer = "".join(true_statements)
        
    print(final_answer)

solve_chromatic_roots_properties()