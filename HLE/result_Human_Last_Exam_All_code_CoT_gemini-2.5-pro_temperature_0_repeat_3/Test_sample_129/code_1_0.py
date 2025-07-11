import itertools

def count_true_expressions():
    """
    Generates all possible expressions of length 5 from the given symbols,
    evaluates them, and counts how many are true.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0

    # Generate all possible strings of length 5
    all_combinations = itertools.product(symbols, repeat=5)

    for combo in all_combinations:
        # Create the expression string, e.g., "T&F|T"
        expr_str = "".join(combo)

        # Convert the expression string into a Python-evaluable format.
        # We add spaces around operators and values to ensure they are treated
        # as distinct tokens by Python's parser.
        py_expr = expr_str.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')
        py_expr = py_expr.replace('!', ' not ')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')

        try:
            # Evaluate the expression in a safe, restricted environment.
            # The 'globals' dictionary only allows 'True' and 'False' to be recognized.
            # The 'locals' dictionary is empty.
            # This prevents the execution of any malicious code.
            result = eval(py_expr, {"True": True, "False": False}, {})
            
            # Check if the expression evaluates to True
            if result is True:
                true_expression_count += 1
        except Exception:
            # If eval() raises any error (e.g., SyntaxError, TypeError),
            # it means the string is not a valid boolean expression.
            # We simply ignore it and move to the next combination.
            pass
            
    return true_expression_count

# Calculate and print the final count.
final_count = count_true_expressions()
print(final_count)