import itertools

def solve_boolean_expressions():
    """
    Finds and counts all true boolean expressions of exactly 5 symbols.
    """
    # The set of allowed symbols for building expressions.
    symbols = ['T', 'F', '!', '&', '|', '(', ')']

    # A list to store the valid, true expressions found.
    true_expressions = []

    # Generate all possible strings of length 5 using the given symbols.
    for p in itertools.product(symbols, repeat=5):
        expr_str = "".join(p)

        # Convert the symbolic expression into a Python-evaluable string.
        # Note: A simple replace works well here.
        py_expr = expr_str.replace('T', 'True')
        py_expr = py_expr.replace('F', 'False')
        py_expr = py_expr.replace('&', 'and')
        py_expr = py_expr.replace('|', 'or')
        # The space after 'not' is crucial to correctly parse expressions like '!T'.
        py_expr = py_expr.replace('!', 'not ')

        try:
            # Use a restricted environment for eval() for security.
            # This ensures only boolean logic is evaluated.
            result = eval(py_expr, {"__builtins__": {}}, {"True": True, "False": False})
            
            # We are only interested in expressions that evaluate to True.
            if result is True:
                true_expressions.append(expr_str)

        except (SyntaxError, TypeError, ValueError, NameError, MemoryError):
            # If eval() fails, the string is not a valid expression.
            # We catch the errors and simply move on to the next string.
            pass
    
    # Print all the true expressions found.
    print("Found the following true expressions:")
    for expr in true_expressions:
        print(expr)

    # Print the final count.
    print(f"\nTotal number of true boolean expressions with exactly 5 symbols: {len(true_expressions)}")

# Run the solver function.
solve_boolean_expressions()