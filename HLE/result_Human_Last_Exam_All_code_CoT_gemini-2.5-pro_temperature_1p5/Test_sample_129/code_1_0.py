import itertools

def count_true_boolean_expressions():
    """
    Generates all possible 5-symbol expressions from the given set,
    evaluates them, and counts how many are true.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0

    # A mapping to translate the symbols into Python's boolean syntax
    translation_map = {
        'T': ' True ',
        'F': ' False ',
        '&': ' and ',
        '|': ' or ',
        '!': ' not ',
        '(': ' ( ',
        ')': ' ) '
    }

    # Generate all possible strings of length 5 with replacement
    for p in itertools.product(symbols, repeat=5):
        # Build the Python expression string from the symbol tuple
        py_expression = "".join(translation_map.get(s, s) for s in p)

        try:
            # First, check if the expression is syntactically valid by compiling it.
            # The 'eval' mode is used for expressions that should return a value.
            code_object = compile(py_expression, '<string>', 'eval')

            # If it compiles without error, evaluate the resulting code object.
            # We specifically check if the result is the boolean True.
            if eval(code_object) is True:
                true_expression_count += 1
        except Exception:
            # Any exception (e.g., SyntaxError, TypeError) indicates an invalid
            # or malformed expression. We ignore these and continue.
            continue

    print(true_expression_count)

if __name__ == '__main__':
    count_true_boolean_expressions()