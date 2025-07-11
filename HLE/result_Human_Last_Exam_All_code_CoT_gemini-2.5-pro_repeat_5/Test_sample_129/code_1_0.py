import itertools

def count_true_expressions():
    """
    Calculates the number of true boolean expressions of a fixed length
    using a given set of symbols.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    length = 5
    true_expression_count = 0

    # Generate all possible strings of the given length
    all_combinations = itertools.product(symbols, repeat=length)

    for combo in all_combinations:
        expression_str = "".join(combo)

        # Translate the expression to a Python-evaluable string
        py_expr = expression_str.replace('T', 'True')
        py_expr = py_expr.replace('F', 'False')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')

        try:
            # Evaluate the expression.
            # We only care about boolean results.
            result = eval(py_expr)
            if isinstance(result, bool) and result is True:
                true_expression_count += 1
        except (SyntaxError, TypeError, NameError):
            # This expression is not syntactically valid, so we ignore it.
            continue
    
    # The final count is the answer to the user's question.
    print(f"Total number of true boolean expressions of length 5 is: {true_expression_count}")

count_true_expressions()