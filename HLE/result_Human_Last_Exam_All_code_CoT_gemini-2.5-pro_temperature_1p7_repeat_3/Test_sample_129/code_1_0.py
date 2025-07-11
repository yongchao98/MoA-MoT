import itertools

def count_true_expressions():
    """
    This function calculates the number of true boolean expressions of length 5
    using the symbols T, F, !, &, |, (, ).
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0

    # Generate all possible strings of length 5. Total permutations = 7^5 = 16807
    all_combinations = itertools.product(symbols, repeat=5)

    for combo in all_combinations:
        expr_str = "".join(combo)

        # Translate the expression into a Python-evaluable format.
        # We add spaces around operators to ensure correct parsing.
        py_expr = expr_str.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')
        # Replace constants after operators to avoid partial replacements like 'False' -> 'F'alse
        py_expr = py_expr.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')

        try:
            # Use eval() to parse and evaluate the expression.
            # The 'is True' check is important to correctly handle non-boolean results,
            # though in this specific grammar, that's unlikely.
            if eval(py_expr) is True:
                true_expression_count += 1
        except (SyntaxError, TypeError, NameError, ValueError, AttributeError, RecursionError):
            # Any string that is not a valid Python expression will raise an error.
            # We catch these errors and simply move on to the next combination.
            pass

    # Print the final count of true expressions.
    print("The final count of true boolean expressions is:")
    print(true_expression_count)

if __name__ == "__main__":
    count_true_expressions()