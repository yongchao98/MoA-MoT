import itertools

def count_true_expressions():
    """
    Generates all possible expressions of length 5 from the given symbols,
    evaluates them, and counts how many are syntactically valid and
    evaluate to True.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    
    # Mapping from expression symbols to Python syntax
    # Note the spaces around the operators to ensure correct parsing
    translation_map = {
        'T': 'True',
        'F': 'False',
        '!': ' not ',
        '&': ' and ',
        '|': ' or ',
        '(': '(',
        ')': ')',
    }

    true_expression_count = 0
    
    # Generate all possible strings of length 5
    for p in itertools.product(symbols, repeat=5):
        expr_str = "".join(p)
        
        # Translate the string into a Python-evaluable expression
        py_expr = expr_str
        py_expr = py_expr.replace('!', translation_map['!'])
        py_expr = py_expr.replace('&', translation_map['&'])
        py_expr = py_expr.replace('|', translation_map['|'])
        py_expr = py_expr.replace('T', translation_map['T'])
        py_expr = py_expr.replace('F', translation_map['F'])
        
        try:
            # Using eval() is safe here because the input strings are generated
            # from a known, safe set of symbols, not from an external user.
            result = eval(py_expr)
            
            # Check if the result is exactly the boolean True
            if result is True:
                true_expression_count += 1
        except (SyntaxError, TypeError, NameError):
            # This expression is not syntactically valid, so we ignore it.
            continue
            
    print(true_expression_count)

if __name__ == "__main__":
    count_true_expressions()