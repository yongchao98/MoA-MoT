import itertools

def find_true_boolean_expressions():
    """
    This function finds, lists, and counts all true boolean expressions of length 5
    that can be formed using the symbols T, F, !, &, |, (, ).
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expressions = []

    # Generate all possible strings of length 5
    all_combinations = itertools.product(symbols, repeat=5)

    for combo in all_combinations:
        expr_str = "".join(combo)

        # Transform to a Python-evaluable expression.
        # Spaces are added to ensure proper tokenization by Python's parser.
        py_expr = expr_str.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')

        # Use a try-except block to safely evaluate the string. This will
        # filter out any syntactically invalid expressions.
        try:
            # A quick pre-check for unbalanced parentheses can speed things up slightly
            # by avoiding some calls to eval that are guaranteed to fail.
            if py_expr.count('(') != py_expr.count(')'):
                continue

            result = eval(py_expr)

            # Check if the result is exactly the boolean value True.
            if result is True:
                true_expressions.append(expr_str)
        except Exception:
            # If eval() fails, it's not a valid expression. We ignore the
            # error and continue to the next combination.
            pass

    # Sort the results for consistent and readable output
    true_expressions.sort()

    # Print the final list of found expressions
    print(f"Found {len(true_expressions)} true boolean expressions of length 5:\n")
    for expr in true_expressions:
        # This loop prints each symbol in the final expressions, fulfilling the prompt.
        print(expr)

if __name__ == '__main__':
    find_true_boolean_expressions()