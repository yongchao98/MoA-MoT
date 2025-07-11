import collections

# Use a dictionary for memoization to cache results for generate_expressions(k)
memo = {}

def generate_expressions(length):
    """
    Generates all syntactically valid boolean expressions of a given length
    using dynamic programming and memoization.
    Returns a set of expression strings.
    """
    # If we have computed this length before, return the cached result.
    if length in memo:
        return memo[length]
    
    # Base cases
    if length <= 0:
        return set()
    if length == 1:
        return {"T", "F"}

    expressions = set()

    # Rule 1: Form '!A' from expressions of length k-1
    if length > 1:
        sub_exprs = generate_expressions(length - 1)
        for s in sub_exprs:
            expressions.add(f"!{s}")

    # Rule 2: Form '(A)' from expressions of length k-2
    if length > 2:
        sub_exprs = generate_expressions(length - 2)
        for s in sub_exprs:
            # We only add parentheses if the sub-expression is not a single literal,
            # to avoid redundant forms like (T). However, the recursive generation
            # might produce complex but valid expressions that need wrapping.
            # A full parser would handle this, but for generation, adding all forms
            # is safer. Duplicates are handled by the set.
            expressions.add(f"({s})")

    # Rule 3: Form 'A op B' from expressions of combined length k-1
    if length > 2:
        # Partition the remaining length (length - 1 for the operator)
        for i in range(1, length - 1):
            j = length - 1 - i
            
            # Get the sets of sub-expressions for the left and right sides
            exprs1 = generate_expressions(i)
            exprs2 = generate_expressions(j)
            
            # Combine them with each operator
            for e1 in exprs1:
                for e2 in exprs2:
                    expressions.add(f"{e1}&{e2}")
                    expressions.add(f"{e1}|{e2}")

    # Cache the result for the current length and return it
    memo[length] = expressions
    return expressions

def evaluate_expression(expr_str):
    """
    Evaluates a boolean expression string by converting it to a Python-evaluable format.
    The precedence is handled correctly by Python's eval: not > and > or.
    """
    try:
        # Replace logic symbols with Python keywords
        py_expr = expr_str.replace('!', ' not ').replace('&', ' and ').replace('|', ' or ')
        py_expr = py_expr.replace('T', ' True ').replace('F', ' False ')
        return eval(py_expr)
    except (SyntaxError, TypeError):
        # This should not be reached if the expression generation is correct
        return None

def solve():
    """
    Main function to find and print true boolean expressions of length 5.
    """
    target_length = 5
    all_expressions = generate_expressions(target_length)
    
    true_expressions = []
    for expr in sorted(list(all_expressions)):
        if evaluate_expression(expr) is True:
            true_expressions.append(expr)

    print(f"Found {len(true_expressions)} true boolean expressions of length {target_length}:")
    
    # The problem asks to "output each number in the final equation"
    # We interpret this as showing each expression that contributes to the final count.
    for expr in true_expressions:
        print(f"{expr}")
        
    print(f"\nTotal true expressions: {len(true_expressions)}")

if __name__ == '__main__':
    solve()
<<<42>>>