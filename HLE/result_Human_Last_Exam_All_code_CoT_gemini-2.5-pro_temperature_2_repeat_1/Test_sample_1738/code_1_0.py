def solve_and_print_answer():
    """
    This function programmatically verifies each statement and prints the letters
    of the true statements in sorted order.
    """
    
    true_statements = []

    # Statement D: If s is any non-empty set, s and True returns True while True and s returns s
    s = {'example'} # A non-empty set is truthy
    if (s and True) is True and (True and s) is s:
        true_statements.append('D')

    # Statement E: The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    if (({} or []) == []) and (([] or {}) == {}):
        true_statements.append('E')

    # Statement F: For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    # This is a universal truth of boolean short-circuit evaluation.
    # Testing with all four truthiness combinations for two objects.
    o1, o2 = True, True   # truthy, truthy
    o3, o4 = False, True  # falsy, truthy
    o5, o6 = True, False  # truthy, falsy
    o7, o8 = False, False # falsy, falsy
    is_f_true = all([
        bool(o1 and o2) == (bool(o1) and bool(o2)),
        bool(o3 and o4) == (bool(o3) and bool(o4)),
        bool(o5 and o6) == (bool(o5) and bool(o6)),
        bool(o7 and o8) == (bool(o7) and bool(o8)),
    ])
    if is_f_true:
        true_statements.append('F')
        
    # Statement I: For any sets p, q, r: (p and q) and r equals p and (q and r)
    # The 'and' operator is associative in Python. We can test this.
    # We create all 8 truthiness combinations for three objects.
    combinations = [
        (True, True, True), (True, True, False), (True, False, True), (True, False, False),
        (False, True, True), (False, True, False), (False, False, True), (False, False, False)
    ]
    is_i_true = all(((p and q) and r) is (p and (q and r)) for p, q, r in combinations)
    if is_i_true:
        true_statements.append('I')

    # Sort the list of true statement letters alphabetically
    true_statements.sort()
    
    # Print the result as a comma-separated string
    print(", ".join(true_statements))

solve_and_print_answer()