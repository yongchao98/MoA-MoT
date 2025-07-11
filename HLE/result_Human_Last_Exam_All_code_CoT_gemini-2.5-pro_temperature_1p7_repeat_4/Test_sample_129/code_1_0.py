import itertools

def count_true_expressions():
    """
    This function calculates the number of true boolean expressions that can be
    written using exactly 5 symbols from the set {T, F, !, &, |, (, )}.

    It works by generating all 7^5 = 16,807 possible strings and evaluating each one.
    """
    count = 0
    symbols = ['T', 'F', '!', '&', '|', '(', ')']

    # The instruction "output each number in the final equation" is interpreted
    # as printing the final result of the count. In this brute-force approach,
    # the count is a single accumulated number.
    final_equation_number = 0

    # Generate all possible expressions of length 5.
    for p in itertools.product(symbols, repeat=5):
        expr_str = "".join(p)

        # Convert the expression into a format Python's eval() can understand.
        # We add spaces around each token to ensure proper parsing.
        # For example, 'T&F' becomes ' True  and  False ', not 'TrueandFalse'.
        py_expr = expr_str.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')
        py_expr = py_expr.replace('(', ' ( ')
        py_expr = py_expr.replace(')', ' ) ')

        # Use a try-except block to filter out syntactically invalid expressions.
        # For example, "T!F&(" would cause a SyntaxError.
        try:
            # An expression must contain a boolean value to be valid.
            # This is a simple check to quickly discard junk like "&&&&&".
            if 'True' not in py_expr and 'False' not in py_expr:
                continue

            # Evaluate the expression. If it's valid and returns True, increment count.
            if eval(py_expr) is True:
                count += 1
        except Exception:
            # If any error occurs during parsing or evaluation, it's not a valid
            # expression we are looking for.
            continue

    final_equation_number = count
    print(final_equation_number)

if __name__ == '__main__':
    count_true_expressions()