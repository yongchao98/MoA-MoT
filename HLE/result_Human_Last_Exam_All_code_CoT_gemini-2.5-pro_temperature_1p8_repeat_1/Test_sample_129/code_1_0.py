import itertools

def count_true_expressions():
    """
    Generates all possible boolean expressions of length 5 from a given
    set of symbols, evaluates them, and counts how many are true.
    """
    # The set of allowed symbols for building expressions.
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    
    # A mapping from our symbols to Python's boolean syntax.
    # Adding spaces around operators prevents syntax errors like 'notTrue'.
    translation_map = {
        'T': ' True ',
        'F': ' False ',
        '&': ' and ',
        '|': ' or ',
        '!': ' not ',
        '(': ' ( ',
        ')': ' ) ',
    }

    true_expression_count = 0
    
    # Use itertools.product to generate all possible expressions of length 5.
    # This creates a generator for all 7^5 = 16,807 combinations.
    all_possible_tuples = itertools.product(symbols, repeat=5)

    for expr_tuple in all_possible_tuples:
        # Translate the tuple of symbols into a Python-evaluable string.
        # e.g., ('(', 'T', '|', 'F', ')') -> ' (  True  or  False  ) '
        py_expr = "".join(translation_map[s] for s in expr_tuple)
        
        try:
            # Use eval() to parse and compute the value of the expression.
            # We wrap this in a try-except block to catch invalid syntax.
            result = eval(py_expr)
            
            # Check if the result is the boolean value True.
            # The isinstance check ensures we only count boolean results,
            # though with these symbols, only booleans are expected.
            if isinstance(result, bool) and result is True:
                true_expression_count += 1
                
        except (SyntaxError, TypeError, NameError, ValueError):
            # These exceptions are expected for malformed expressions.
            # e.g., "T & & T" would raise a SyntaxError.
            # We simply ignore invalid expressions and continue.
            pass
            
    # Print the final count, which is the answer to the problem.
    print(true_expression_count)

# Execute the function to find the answer.
count_true_expressions()