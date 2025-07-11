import itertools

def solve_boolean_expression_count():
    """
    Calculates and prints the number of true boolean expressions
    that can be written using exactly 5 symbols from the given set.
    """
    # The set of allowed symbols for building expressions.
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0

    # A mapping from our custom symbols to Python's boolean syntax.
    symbol_to_python_map = {
        'T': 'True',
        'F': 'False',
        '!': 'not',
        '&': 'and',
        '|': 'or',
        '(': '(',
        ')': ')'
    }

    # Generate all possible combinations of 5 symbols, with repetition allowed.
    all_combinations = itertools.product(symbols, repeat=5)

    # Iterate through each generated combination.
    for combo in all_combinations:
        # Optimization: Perform quick checks to filter out many obviously invalid expressions.
        # An expression cannot start with a binary operator or a closing parenthesis.
        if combo[0] in ('&', '|', ')'):
            continue
        # An expression cannot end with any operator or an opening parenthesis.
        if combo[-1] in ('&', '|', '!', '('):
            continue
        # The number of opening and closing parentheses must be equal for a valid expression.
        if combo.count('(') != combo.count(')'):
            continue

        # Convert the combination tuple into a list of Python-equivalent tokens.
        # e.g., ('!', '!', 'T', '&', 'T') becomes ['not', 'not', 'True', 'and', 'True'].
        python_tokens = [symbol_to_python_map[s] for s in combo]

        # Join the tokens with spaces to form a valid Python expression string.
        # e.g., 'not not True and True'
        python_expression_string = ' '.join(python_tokens)

        try:
            # Use eval() to parse and evaluate the expression string.
            # A try-except block is used to catch any errors from malformed expressions.
            result = eval(python_expression_string)

            # We are only interested in expressions that evaluate to the boolean True.
            if result is True:
                true_expression_count += 1
        except Exception:
            # If any error occurs during evaluation, it's not a valid expression.
            # Ignore it and continue to the next combination.
            continue

    # As requested by the prompt "How many...", we print the final count.
    # The instruction "output each number in the final equation" is interpreted
    # as printing the single final number of the count.
    print(true_expression_count)

solve_boolean_expression_count()