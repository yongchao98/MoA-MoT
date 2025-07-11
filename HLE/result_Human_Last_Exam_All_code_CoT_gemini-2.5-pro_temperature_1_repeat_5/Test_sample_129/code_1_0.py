import itertools

def find_true_boolean_expressions():
    """
    This function generates all possible expressions of length 5 from a given
    set of symbols, evaluates them, and prints those that are syntactically
    valid and evaluate to True.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    length = 5
    
    true_expressions = []
    
    # Generate all possible combinations of symbols of the given length
    all_combinations = itertools.product(symbols, repeat=length)
    
    for combo in all_combinations:
        # Join the characters to form the expression string
        expr_str = "".join(combo)
        
        # Convert the custom expression syntax to Python's boolean syntax
        # We add spaces around operators to ensure correct parsing by eval()
        py_expr = expr_str.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')
        py_expr = py_expr.replace('T', 'True')
        py_expr = py_expr.replace('F', 'False')
        
        try:
            # Safely evaluate the expression in a restricted environment.
            # Only 'True' and 'False' are defined in the context.
            result = eval(py_expr, {"__builtins__": {}}, {"True": True, "False": False})
            
            # Check if the expression evaluates to the boolean value True
            if result is True:
                true_expressions.append(expr_str)
        except Exception:
            # Any string that causes an error (e.g., SyntaxError) is not a
            # valid expression, so we ignore it and continue.
            pass
            
    # As requested, output each expression found
    print("The true boolean expressions of length 5 are:")
    for expr in true_expressions:
        print(expr)
        
    print("\nTotal number of true expressions found:")
    print(len(true_expressions))

# Run the function to find and print the expressions and the final count
find_true_boolean_expressions()
<<<81>>>