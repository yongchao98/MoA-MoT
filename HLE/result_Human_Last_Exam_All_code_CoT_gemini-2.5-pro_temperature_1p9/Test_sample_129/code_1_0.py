import itertools

def is_valid_expression(expr_str: str) -> bool:
    """
    Checks if a given string is a valid boolean expression.
    - Translates the custom symbols to Python's boolean operators.
    - Uses eval() in a restricted context to parse and evaluate.
    - A valid expression must evaluate to a boolean without raising errors.
    """
    # A quick and necessary check for balanced parentheses.
    if expr_str.count('(') != expr_str.count(')'):
        return False
    
    # Replace symbols with Python's keywords. Note the spaces around keywords
    # to prevent incorrect concatenations (e.g., 'TF' becoming 'TrueFalse').
    py_expr = expr_str.replace('T', ' True ')
    py_expr = py_expr.replace('F', ' False ')
    py_expr = py_expr.replace('&', ' and ')
    py_expr = py_expr.replace('|', ' or ')
    py_expr = py_expr.replace('!', ' not ')
    
    try:
        # Use a restricted global and local namespace for eval safety.
        # __builtins__ is empty to prevent access to built-in functions.
        # Only 'True' and 'False' are defined as names.
        result = eval(py_expr, {'__builtins__': {}}, {'True': True, 'False': False})
        
        # Ensure the final result is a boolean value.
        if isinstance(result, bool):
            return True
        return False
    except Exception:
        # Any exception during parsing or evaluation means it's not a valid expression.
        return False

def solve_boolean_expressions():
    """
    Finds and counts all valid boolean expressions of exactly 5 symbols.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    
    # A dictionary to hold the counts for each category, based on the starting symbol.
    counts_by_start = {
        'T': 0,
        'F': 0,
        '!': 0,
        '(': 0,
    }

    # Generate all possible strings of length 5.
    for p in itertools.product(symbols, repeat=5):
        expression = "".join(p)
        
        # Check if the generated expression is valid.
        if is_valid_expression(expression):
            first_char = expression[0]
            if first_char in counts_by_start:
                counts_by_start[first_char] += 1
    
    # Prepare the components for the final "equation" output.
    equation_parts = []
    for symbol in ['T', 'F', '!', '(']:
        count = counts_by_start[symbol]
        if count > 0:
            equation_parts.append(str(count))
            
    total_count = sum(counts_by_start.values())
    
    # Print the breakdown as an equation.
    equation_str = " + ".join(equation_parts)
    print(f"Found {total_count} valid expressions.")
    print(f"Breakdown by starting symbol: {equation_str} = {total_count}")


if __name__ == "__main__":
    solve_boolean_expressions()