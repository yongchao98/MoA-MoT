import itertools

def solve_boolean_expressions():
    """
    Finds and counts all true boolean expressions of a fixed length
    using a given set of symbols.
    """
    # The set of allowed symbols for the expressions.
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    
    # We will generate all possible strings of length 5.
    expression_length = 5
    all_possible_expressions = itertools.product(symbols, repeat=expression_length)
    
    true_expression_count = 0
    
    # Iterate through each generated expression string.
    for expr_tuple in all_possible_expressions:
        expression_str = "".join(expr_tuple)
        
        # Convert the expression into a format Python can evaluate.
        # We add spaces to ensure tokens are separated, e.g., 'T&T' -> ' True  and  True '.
        py_expr = expression_str.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')

        try:
            # Evaluate the expression in a controlled, safe environment.
            # Only 'True' and 'False' are defined, and no built-ins are available.
            result = eval(py_expr, {'__builtins__': None}, {'True': True, 'False': False})
            
            # Check if the expression is syntactically valid and evaluates to True.
            if result is True:
                true_expression_count += 1
                
        except Exception:
            # Any exception (e.g., SyntaxError, TypeError) means the expression is invalid.
            # We simply ignore it and move to the next one.
            continue
            
    # Print the final results.
    print(f"Found {true_expression_count} true boolean expressions of length 5.")
    
    # To satisfy the instruction "output each number in the final equation",
    # we represent the total count as a sum of 1s.
    print("\nThis result can be shown with the following equation, where each '1' represents one valid true expression:")
    if true_expression_count > 0:
        # To keep the output tidy, we don't print thousands of '1's.
        max_ones_to_show = 10
        if true_expression_count <= max_ones_to_show:
            sum_of_ones = " + ".join(["1"] * true_expression_count)
        else:
            sum_of_ones = " + ".join(["1"] * max_ones_to_show) + " + ..."
            
        print(f"{sum_of_ones} = {true_expression_count}")
    else:
        print("0 = 0")

if __name__ == '__main__':
    solve_boolean_expressions()