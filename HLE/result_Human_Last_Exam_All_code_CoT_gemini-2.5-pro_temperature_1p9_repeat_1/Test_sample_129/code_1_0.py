import itertools

def find_true_boolean_expressions():
    """
    This function generates all possible 5-symbol boolean expressions from a given
    set of symbols, evaluates them, and prints those that result in True.
    """
    # The set of allowed symbols
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    
    true_expressions = []
    
    # 1. Generate all possible combinations of 5 symbols
    total_combinations = itertools.product(symbols, repeat=5)
    
    for combo in total_combinations:
        # A quick check for unbalanced parentheses, which are always invalid
        if combo.count('(') != combo.count(')'):
            continue

        original_expr = "".join(combo)
        
        # 2. Translate the expression into a Python-evaluable string.
        # Spacing around operators is crucial for Python's parser.
        py_expr = original_expr.replace('!', 'not ') \
                               .replace('&', ' and ') \
                               .replace('|', ' or ') \
                               .replace('T', 'True') \
                               .replace('F', 'False')
        
        # 3. Use a try-except block to safely evaluate the string
        try:
            # The 'eval' function will parse and compute the result.
            # It will raise an error for malformed expressions (e.g., "T&&F").
            # We explicitly check for `is True` to ensure the result is a boolean True.
            if eval(py_expr) is True:
                true_expressions.append(original_expr)
        except Exception:
            # This catches SyntaxError, TypeError, etc. for invalid expressions.
            # We simply ignore them and move to the next combination.
            pass

    # 4. Print all the valid true expressions found, sorted for consistency.
    print("Found the following true boolean expressions of length 5:")
    for expr in sorted(true_expressions):
        print(expr)
        
    print(f"\nTotal number of true expressions: {len(true_expressions)}")

# Execute the function
find_true_boolean_expressions()