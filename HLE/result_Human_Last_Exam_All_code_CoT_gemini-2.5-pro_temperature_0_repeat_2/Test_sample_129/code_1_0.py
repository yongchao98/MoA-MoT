import itertools

def find_true_boolean_expressions():
    """
    This function finds all true boolean expressions of exactly 5 symbols
    from the set {T, F, !, &, |, (, )}, prints them as numeric equations,
    and returns the total count.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0
    
    print("Finding all true boolean expressions of length 5...")
    print("The expressions are shown below in numeric form (T=1, F=0):")
    print("-" * 30)

    # 1. Generate all possible strings of length 5
    all_combinations = itertools.product(symbols, repeat=5)
    
    for combo in all_combinations:
        expression_str = "".join(combo)
        
        # 2. Convert to a Python-evaluatable string
        # Add spaces around operators to ensure correct parsing by eval()
        py_expr = expression_str.replace('T', ' True ')
        py_expr = py_expr.replace('F', ' False ')
        py_expr = py_expr.replace('&', ' and ')
        py_expr = py_expr.replace('|', ' or ')
        py_expr = py_expr.replace('!', ' not ')

        try:
            # 3. Evaluate the expression and check if it is True
            result = eval(py_expr)
            if result is True:
                true_expression_count += 1
                
                # 4. Format and print the numeric equation as requested
                # This fulfills the instruction "output each number in the final equation!"
                numeric_expr = expression_str.replace('T', '1').replace('F', '0')
                print(f"{numeric_expr} = 1")

        except (SyntaxError, TypeError, NameError, ValueError):
            # If eval fails, it's not a valid boolean expression, so we ignore it.
            pass
            
    # 5. Print the final total count
    print("-" * 30)
    print(f"The total number of true boolean expressions is: {true_expression_count}")

if __name__ == '__main__':
    find_true_boolean_expressions()