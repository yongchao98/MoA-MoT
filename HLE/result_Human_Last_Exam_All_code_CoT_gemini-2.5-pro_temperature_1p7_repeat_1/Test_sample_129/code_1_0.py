import itertools

def count_true_expressions():
    """
    Generates all valid boolean expressions of exactly 5 symbols,
    evaluates them, and counts how many are true.
    """
    true_expressions = set()
    operands = ['T', 'F']
    binary_ops = ['&', '|']

    # Helper for converting to a Python-evaluable string
    def to_py_expr(s):
        return s.replace('T', 'True').replace('F', 'False').replace('&', 'and').replace('|', 'or').replace('!', 'not ')

    # Pattern 1: A op1 B op2 C
    for a, b, c in itertools.product(operands, repeat=3):
        for op1, op2 in itertools.product(binary_ops, repeat=2):
            expr_str = f"{a}{op1}{b}{op2}{c}"
            if eval(to_py_expr(expr_str)):
                true_expressions.add(expr_str)

    # Pattern 2: !!A op B
    for a, b in itertools.product(operands, repeat=2):
        for op in binary_ops:
            expr_str = f"!!{a}{op}{b}"
            if eval(to_py_expr(expr_str)):
                true_expressions.add(expr_str)
    
    # Pattern 3: A op !!B
    for a, b in itertools.product(operands, repeat=2):
        for op in binary_ops:
            expr_str = f"{a}{op}!!{b}"
            if eval(to_py_expr(expr_str)):
                true_expressions.add(expr_str)

    # Pattern 4: !A op !B
    for a, b in itertools.product(operands, repeat=2):
        for op in binary_ops:
            expr_str = f"!{a}{op}!{b}"
            if eval(to_py_expr(expr_str)):
                true_expressions.add(expr_str)

    # Pattern 5: !!!!A
    for a in operands:
        expr_str = f"!!!!{a}"
        if eval(to_py_expr(expr_str)):
            true_expressions.add(expr_str)
    
    # Pattern 6: (A op B)
    for a, b in itertools.product(operands, repeat=2):
        for op in binary_ops:
            expr_str = f"({a}{op}{b})"
            if eval(to_py_expr(expr_str)):
                true_expressions.add(expr_str)

    # Pattern 7: (!!A)
    for a in operands:
        expr_str = f"(!!{a})"
        if eval(to_py_expr(expr_str)):
            true_expressions.add(expr_str)
    
    # Pattern 8: ((A))
    for a in operands:
        expr_str = f"(({a}))"
        if eval(to_py_expr(expr_str)):
            true_expressions.add(expr_str)

    # Pattern 9: !(!A)
    for a in operands:
        expr_str = f"!(!{a})"
        if eval(to_py_expr(expr_str)):
            true_expressions.add(expr_str)
    
    # Pattern 10: (A) op B
    for a, b in itertools.product(operands, repeat=2):
        for op in binary_ops:
            expr_str = f"({a}){op}{b}"
            if eval(to_py_expr(expr_str)):
                true_expressions.add(expr_str)
    
    # Pattern 11: A op (B)
    for a, b in itertools.product(operands, repeat=2):
        for op in binary_ops:
            expr_str = f"{a}{op}({b})"
            if eval(to_py_expr(expr_str)):
                true_expressions.add(expr_str)
    
    # Print all the found true expressions
    print("Found the following true boolean expressions:")
    for expr in sorted(list(true_expressions)):
        print(expr)
    
    # Print the final count
    print(f"\nTotal number of unique true expressions is: {len(true_expressions)}")

if __name__ == '__main__':
    count_true_expressions()
<<<46>>>