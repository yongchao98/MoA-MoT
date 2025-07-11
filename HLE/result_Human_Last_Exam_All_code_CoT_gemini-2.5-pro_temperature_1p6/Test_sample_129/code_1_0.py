import itertools

def count_true_expressions():
    """
    Calculates the number of true boolean expressions that can be written
    using exactly 5 symbols from the set {T, F, !, &, |, (, )}.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0

    # This mapping translates the given symbols into valid Python syntax.
    # Operators are padded with spaces to ensure they are treated as separate tokens.
    # Operands (T, F) and parentheses do not need padding.
    translation_map = {
        'T': 'True',
        'F': 'False',
        '!': ' not ',
        '&': ' and ',
        '|': ' or ',
        '(': '(',
        ')': ')'
    }

    # Generate all possible permutations of length 5 with repetition.
    for p in itertools.product(symbols, repeat=5):
        expression_str = "".join(p)

        # Build the Python-evaluatable string using the translation map.
        py_expr_tokens = [translation_map[token] for token in expression_str]
        py_expr = "".join(py_expr_tokens)

        try:
            # eval() parses and evaluates the string. This is a powerful but
            # potentially unsafe function. It's used here because we are generating
            # the input strings ourselves from a controlled set of symbols.
            result = eval(py_expr)

            # We are only interested in expressions that are explicitly True.
            if result is True:
                true_expression_count += 1
        except (SyntaxError, TypeError, NameError, ValueError, AttributeError):
            # If eval() fails, the string is not a valid expression (e.g., "T&&()"),
            # so we ignore it and continue to the next possibility.
            pass

    # The problem asks for the total count. The "equation" is the final count itself.
    print("Total number of true expressions of length 5 is:", true_expression_count)

if __name__ == "__main__":
    count_true_expressions()