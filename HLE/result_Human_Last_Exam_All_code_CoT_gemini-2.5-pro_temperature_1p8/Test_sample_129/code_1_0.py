import itertools

def count_true_expressions():
    """
    This function calculates the number of true boolean expressions of length 5
    using the symbols T, F, !, &, |, (, ). It works by generating all possible
    strings and using Python's eval() to check for syntactic validity and truth value.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    
    # Using a set to store unique expressions that evaluate to True
    true_expressions = set()

    # The total number of combinations to check is len(symbols)^5
    total_combinations = len(symbols) ** 5

    # Generate all possible strings of length 5
    for p in itertools.product(symbols, repeat=5):
        expr_str = "".join(p)

        # Convert the expression to a Python-evaluable format.
        # Spaces are added around keywords to ensure proper tokenization.
        py_expr_str = expr_str.replace('T', ' True ')
        py_expr_str = py_expr_str.replace('F', ' False ')
        py_expr_str = py_expr_str.replace('&', ' and ')
        py_expr_str = py_expr_str.replace('|', ' or ')
        py_expr_str = py_expr_str.replace('!', ' not ')

        try:
            # The 'eval' function parses and evaluates the expression.
            # It naturally handles operator precedence (not > and > or) and syntax.
            # Malformed expressions (e.g., "T & & F", "( T ) F") will raise an exception.
            result = eval(py_expr_str)

            if result is True:
                # Add the original valid expression to our set of solutions.
                true_expressions.add(expr_str)
        except Exception:
            # This catches SyntaxError, TypeError, etc., for invalid expressions.
            # We simply ignore them and continue to the next combination.
            pass
            
    # The final answer is the number of unique true expressions found.
    # The problem asks for the count, which is a single number.
    final_count = len(true_expressions)
    
    # As per the instructions, we output the "final equation" numbers.
    # Here, the process is an exhaustive search, not a simple equation.
    # We will output the final count as the result of our work.
    print(final_count)

if __name__ == '__main__':
    count_true_expressions()