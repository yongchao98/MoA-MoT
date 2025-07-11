def solve_chromatic_roots_properties():
    """
    Analyzes statements about chromatic roots and identifies the true ones.
    
    The truth value of each statement is determined based on well-known theorems in graph theory.
    A: True - Real orbital chromatic roots are bounded by the greatest real chromatic root (Cameron & Kayibi, 2007).
    B: True - Chromatic roots can be complex. e.g., for C4, P(k) = k(k-1)(k^2-3k+3) has complex roots.
    C: False - Tutte showed P(G, k) is non-zero for all k < 0.
    D: True - Tutte showed that non-integer roots exist, e.g., φ^2 ≈ 2.618.
    E: False - Tutte showed P(G, k) is non-zero for all k in the interval (0, 1).
    """
    
    # Store the truth values for each statement
    statement_truth_values = {
        'A': True,
        'B': True,
        'C': False,
        'D': True,
        'E': False
    }
    
    true_statements = []
    for statement, is_true in statement_truth_values.items():
        if is_true:
            true_statements.append(statement)
            
    # Sort the letters alphabetically
    sorted_true_statements = sorted(true_statements)
    
    # Join the letters to form the final answer string
    final_answer = "".join(sorted_true_statements)
    
    print(final_answer)

solve_chromatic_roots_properties()