import itertools

def count_true_expressions():
    """
    Calculates the number of true boolean expressions of length 5 using the symbols
    T, F, !, &, |, (, ).
    The code iterates through all 7^5 = 16,807 possible strings, transforms them
    into a Python-evaluable format, and safely evaluates them.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0

    # Generate all possible strings of length 5.
    for p in itertools.product(symbols, repeat=5):
        expression_str = "".join(p)

        # Transform the expression string into a Python-evaluable format.
        # Adding spaces around all symbols ensures that they are treated as separate tokens.
        # For example, '!T' becomes ' not  True ', preventing 'notTrue' which would be a NameError.
        py_expr_str = expression_str.replace('T', ' True ') \
                                    .replace('F', ' False ') \
                                    .replace('&', ' and ') \
                                    .replace('|', ' or ') \
                                    .replace('!', ' not ') \
                                    .replace('(', ' ( ') \
                                    .replace(')', ' ) ')

        try:
            # The context `{"__builtins__": {}}` is used for safety, to prevent the expression
            # from accessing any built-in functions. The second dictionary `{}` is for local variables.
            result = eval(py_expr_str, {"__builtins__": {}}, {})

            # The problem asks for *boolean* expressions. We must ensure the evaluated
            # result is strictly the boolean value True. This check prevents counting
            # "truthy" values of other types (e.g., non-empty tuples or non-zero numbers).
            if result is True:
                true_expression_count += 1
        except Exception:
            # Any exception (SyntaxError, NameError, TypeError, etc.) means the string
            # is not a syntactically correct and evaluatable boolean expression.
            # We simply ignore such strings and continue to the next one.
            continue
            
    # The phrase "output each number in the final equation" from the prompt is
    # interpreted here as providing the final calculated count.
    print(true_expression_count)

if __name__ == "__main__":
    count_true_expressions()