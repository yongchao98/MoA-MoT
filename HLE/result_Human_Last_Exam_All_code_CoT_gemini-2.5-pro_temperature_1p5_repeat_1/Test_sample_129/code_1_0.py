import itertools

def solve_boolean_expressions():
    """
    Calculates the number of true boolean expressions of exactly length 5
    using the symbols T, F, !, &, |, (, ).
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0
    
    # Generate all possible strings of length 5
    all_combinations = itertools.product(symbols, repeat=5)
    
    for combo in all_combinations:
        expression_str = "".join(combo)
        
        # Translate the string into a valid Python expression
        # The spaces are important for Python's syntax (e.g., 'True and False' not 'TrueandFalse')
        py_expr = expression_str.replace('T', 'True') \
                                .replace('F', 'False') \
                                .replace('&', ' and ') \
                                .replace('|', ' or ') \
                                .replace('!', ' not ')

        try:
            # Use eval() to compute the value of the expression.
            # We check if the result is specifically the boolean True.
            # The try-except block handles all invalid syntax combinations gracefully.
            result = eval(py_expr, {}, {})
            if result is True:
                true_expression_count += 1
        except (SyntaxError, TypeError, NameError, ValueError):
            # If the expression is not valid Python syntax, ignore it.
            continue
            
    # The final count is the answer.
    # The prompt requests to "output each number in the final equation".
    # Since there is no equation, we will print the single resulting number.
    print(f"The total number of true boolean expressions is: {true_expression_count}")

solve_boolean_expressions()