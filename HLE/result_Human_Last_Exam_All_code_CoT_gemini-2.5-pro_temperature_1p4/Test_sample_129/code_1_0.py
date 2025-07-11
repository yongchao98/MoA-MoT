import itertools

def count_true_expressions():
    """
    Generates all possible boolean expressions of a fixed length and counts how many
    evaluate to True.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    length = 5
    true_count = 0

    # A dictionary to map the given symbols to Python's boolean syntax.
    # Note the spaces are crucial for correct parsing (e.g., 'T&F' becomes 'True and False').
    py_syntax_map = {
        'T': 'True',
        'F': 'False',
        '!': 'not ',
        '&': ' and ',
        '|': ' or ',
        '(': '(',
        ')': ')',
    }

    # Generate all possible combinations of 5 symbols with repetitions.
    all_combinations = itertools.product(symbols, repeat=length)

    true_expressions_found = []

    for combo in all_combinations:
        # Create the original expression string (e.g., "(T|F)")
        original_expr = "".join(combo)
        
        # Create the Python-evaluable expression string (e.g., "(True or False)")
        python_expr = "".join(py_syntax_map[symbol] for symbol in combo)

        try:
            # Safely evaluate the expression string.
            # We check for `is True` to be precise and avoid "truthy" values.
            if eval(python_expr) is True:
                true_expressions_found.append(original_expr)
                true_count += 1
        except (SyntaxError, TypeError, NameError):
            # Ignore any expressions that are not syntactically valid.
            # - SyntaxError for bad structure (e.g., "T & & F T").
            # - TypeError for misusing operators (e.g., "not ()").
            # - NameError for things like "True True".
            pass

    # As requested, output each expression that results in True.
    for expr in true_expressions_found:
        print(expr)
    
    # Finally, print the total count.
    print(f"\nTotal number of true expressions: {true_count}")

if __name__ == "__main__":
    count_true_expressions()