import itertools

def count_true_boolean_expressions():
    """
    Calculates the number of true boolean expressions of length 5.

    This function generates all possible 5-symbol strings from the set
    {T, F, !, &, |, (, )}, checks each for syntactic validity as a
    boolean expression, and evaluates the valid ones. It then prints a
    summary of the counts.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    length = 5

    true_count = 0
    false_count = 0

    # Generate all possible strings of the given length
    all_combinations = itertools.product(symbols, repeat=length)

    for combo in all_combinations:
        expr_str = "".join(combo)

        # Convert to a Python-evaluatable string
        py_expr = expr_str.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        # 'not' is a prefix operator, space after is important
        py_expr = py_expr.replace('!', ' not ')

        try:
            # Evaluate the expression
            result = eval(py_expr)

            # Check if the result is a boolean
            if isinstance(result, bool):
                if result is True:
                    true_count += 1
                else:
                    false_count += 1
        except (SyntaxError, TypeError, NameError):
            # Ignore invalid expressions
            pass

    valid_count = true_count + false_count
    
    # Print the "final equation" as requested
    print(f"Found {true_count} true expressions and {false_count} false expressions.")
    print(f"This makes a total of {valid_count} valid boolean expressions.")
    print(f"The equation is: {true_count} + {false_count} = {valid_count}")
    print(f"\nThe number of true boolean expressions is: {true_count}")


if __name__ == "__main__":
    count_true_boolean_expressions()

<<<450>>>