import itertools

def count_true_expressions():
    """
    Counts the number of true boolean expressions of length 5.
    """
    # The set of allowed symbols for building the expression.
    symbols = 'TF!&|()'
    true_expression_count = 0

    # Generate all possible strings of length 5 using the given symbols.
    all_possible_strings = itertools.product(symbols, repeat=5)

    for p in all_possible_strings:
        expression_string = "".join(p)

        # Convert the custom expression syntax to Python's boolean syntax.
        # We add spaces around 'and'/'or' and after 'not' to ensure proper parsing
        # of concatenated expressions like 'T&F' or '!T'.
        python_expression = expression_string.replace('!', ' not ')
        python_expression = python_expression.replace('&', ' and ')
        python_expression = python_expression.replace('|', ' or ')
        python_expression = python_expression.replace('T', 'True')
        python_expression = python_expression.replace('F', 'False')

        try:
            # Evaluate the string as a Python expression.
            # A `try-except` block handles the many syntactically invalid strings.
            result = eval(python_expression)

            # Check if the expression evaluates to True. We use 'is True' for a strict check.
            if result is True:
                true_expression_count += 1
        except Exception:
            # Catch any exceptions from eval() (e.g., SyntaxError, TypeError)
            # and simply continue to the next string.
            continue

    # Print the final result.
    print(f"The total number of true boolean expressions of length 5 is: {true_expression_count}")

count_true_expressions()