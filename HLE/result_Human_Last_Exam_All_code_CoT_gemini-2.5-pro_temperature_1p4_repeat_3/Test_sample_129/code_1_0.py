import itertools

def count_boolean_expressions():
    """
    Counts the number of valid boolean expressions of length 5 using a given set of symbols.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    length = 5
    
    valid_expressions = []
    
    # 1. Generate all possible strings of length 5.
    all_combinations = itertools.product(symbols, repeat=length)
    
    for combo in all_combinations:
        expr_str = "".join(combo)
        
        # 2. Transform the string into a Python-evaluable format.
        # Adding spaces around operators and literals prevents invalid concatenations.
        py_expr = expr_str.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')
        
        # 3. Use eval() in a try-except block to check for syntactic validity.
        try:
            # An expression is valid if it can be evaluated without error
            # and results in a boolean value.
            result = eval(py_expr)
            if isinstance(result, bool):
                valid_expressions.append(expr_str)
        except (SyntaxError, TypeError, NameError):
            # These exceptions indicate an invalid expression format.
            pass
            
    # 4. Categorize the valid expressions and print the results.
    expr_with_parens = []
    expr_without_parens = []
    
    for expr in valid_expressions:
        if '(' in expr or ')' in expr:
            expr_with_parens.append(expr)
        else:
            expr_without_parens.append(expr)
            
    count_with_parens = len(expr_with_parens)
    count_without_parens = len(expr_without_parens)
    total_count = len(valid_expressions)
    
    print("Analysis of valid boolean expressions of length 5:")
    print(f"Found {count_without_parens} expressions without parentheses.")
    print(f"Found {count_with_parens} expressions with parentheses.")
    print("\nFinal count expressed as an equation:")
    
    # Output each number in the final equation as requested.
    print(f"{count_without_parens} + {count_with_parens} = {total_count}")

# Run the function to get the answer.
count_boolean_expressions()
<<<98>>>