import itertools

def count_true_boolean_expressions():
    """
    This script finds the number of true boolean expressions of length 5.

    It generates all possible 5-character strings from the given symbols,
    evaluates them as boolean expressions, and counts how many are true.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    length = 5
    true_expression_count = 0

    # 1. Generate all possible strings of the given length (7^5 = 16,807).
    all_possible_strings = itertools.product(symbols, repeat=length)

    for s_tuple in all_possible_strings:
        expr_str = "".join(s_tuple)

        # 2. Transform the string into a Python-evaluable format.
        # Spaces are added to ensure symbols are parsed as distinct tokens.
        eval_str = expr_str.replace('T', ' True ')
        eval_str = eval_str.replace('F', ' False ')
        eval_str = eval_str.replace('&', ' and ')
        eval_str = eval_str.replace('|', ' or ')
        eval_str = eval_str.replace('!', ' not ')

        try:
            # 3. Evaluate the transformed string safely.
            # We check if the result is specifically the boolean value True.
            if eval(eval_str) is True:
                true_expression_count += 1
        except Exception:
            # This catches any syntactically invalid expressions.
            # We simply ignore them and continue to the next string.
            pass

    # 4. Print the final count.
    # The prompt requests to "output each number in the final equation".
    # For this counting problem, this is interpreted as printing the final total.
    print(true_expression_count)

# Execute the function to find and print the answer.
count_true_boolean_expressions()