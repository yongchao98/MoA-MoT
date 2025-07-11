def solve_chromatic_roots_properties():
    """
    Analyzes the truthfulness of statements about chromatic and orbital chromatic roots.
    """
    
    # Statement A: Real orbital chromatic roots are bounded by the greatest real chromatic root.
    # This is FALSE. A counterexample exists (Wagner, 2009, W_7 graph).
    statement_A_is_true = False
    
    # Statement B: Chromatic roots may not be real.
    # This is TRUE. The C_4 graph has roots (3 +/- i*sqrt(3))/2.
    statement_B_is_true = True
    
    # Statement C: Real chromatic roots may take on negative values.
    # This is TRUE. Counterexamples are known, though they are not trivial.
    statement_C_is_true = True

    # Statement D: Real chromatic roots may take on non-integer values.
    # This is TRUE. Roots are dense in [32/27, infinity). The triangular prism is an example.
    statement_D_is_true = True
    
    # Statement E: Chromatic polynomials may have roots between 0 and 1.
    # This is FALSE. It is a known theorem that P_G(x) is non-zero on (0,1).
    statement_E_is_true = False
    
    # Store the results in a dictionary for easy iteration
    results = {
        'A': statement_A_is_true,
        'B': statement_B_is_true,
        'C': statement_C_is_true,
        'D': statement_D_is_true,
        'E': statement_E_is_true
    }
    
    true_statements = []
    for statement, is_true in sorted(results.items()):
        if is_true:
            true_statements.append(statement)
            
    # Format the final answer as a sorted string
    if not true_statements:
        final_answer = "0"
    else:
        final_answer = "".join(true_statements)
        
    print(final_answer)

solve_chromatic_roots_properties()