import itertools

def evaluate_statements():
    """
    Evaluates each statement from the problem description to determine if it's true or false.
    """
    true_statements = []
    
    # --- Statement A ---
    print("\n--- Evaluating A ---")
    print("Statement: The expression c[0] and c[1] always returns the same value as c[1] and c[0]")
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    expr1 = c[0] and c[1]
    expr2 = c[1] and c[0]
    print(f"Given a = {a}, b = {b}, c = (a, b)")
    print(f"c[0] is truthy, so 'c[0] and c[1]' evaluates to c[1]: {expr1}")
    print(f"c[1] is truthy, so 'c[1] and c[0]' evaluates to c[0]: {expr2}")
    is_true = (expr1 == expr2)
    print(f"Result: The values are not the same. Statement A is False.")
    if is_true: true_statements.append('A')

    # --- Statement B ---
    print("\n--- Evaluating B ---")
    print("Statement: For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2")
    s1 = {1, 2}
    s2 = {2, 3}
    print(f"Counterexample: s1 = {s1}, s2 = {s2}")
    print("s1 and s2 (s1 is truthy) -> s2")
    print("s2 and s1 (s2 is truthy) -> s1")
    result = (s1 and s2) | (s2 and s1)
    print(f"The expression becomes s2 | s1, which is the union: {result}")
    is_true = (result == s2)
    print(f"Result: {result} is not equal to s2 ({s2}). Statement B is False.")
    if is_true: true_statements.append('B')
    
    # --- Statement C ---
    print("\n--- Evaluating C ---")
    print("Statement: The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples")
    d_counter = ((), (3, 4))
    print("The phrase 'regardless of the values' implies it should hold for any tuples, including empty ones.")
    print(f"Counterexample: d = {d_counter}")
    result = d_counter[0] or d_counter[1]
    print(f"d[0] is an empty tuple, which is falsy. So, 'd[0] or d[1]' evaluates to d[1]: {result}")
    is_true = (result == (1, 2))
    print(f"Result: The result is not (1, 2). Statement C is False.")
    if is_true: true_statements.append('C')

    # --- Statement D ---
    print("\n--- Evaluating D ---")
    print("Statement: If s is any non-empty set, s and True returns True while True and s returns s")
    s = {"any", "non-empty", "set"}
    print(f"Test with s = {s}. s is truthy.")
    res1 = s and True
    res2 = True and s
    print(f"'s and True' evaluates to the second operand, True. Result: {res1}")
    print(f"'True and s' evaluates to the second operand, s. Result: {res2}")
    is_true = (res1 is True) and (res2 is s)
    print("Result: Both parts of the statement hold. Statement D is True.")
    if is_true: true_statements.append('D')

    # --- Statement E ---
    print("\n--- Evaluating E ---")
    print("Statement: The expression {} or [] evaluates to [] while [] or {} evaluates to {}")
    res1 = {} or []
    res2 = [] or {}
    print(f"{} is falsy, so '{} or []' evaluates to the second operand: {res1}")
    print(f"[] is falsy, so '[] or {{}}' evaluates to the second operand: {res2}")
    is_true = (res1 == []) and (res2 == {})
    print("Result: Both parts of the statement hold. Statement E is True.")
    if is_true: true_statements.append('E')
    
    # --- Statement F ---
    print("\n--- Evaluating F ---")
    print("Statement: For any tuple t containing two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])")
    s_truthy = {1}
    s_falsy = set()
    cases = [(s_truthy, s_truthy), (s_truthy, s_falsy), (s_falsy, s_truthy), (s_falsy, s_falsy)]
    all_cases_hold = True
    for s1, s2 in cases:
        t = (s1, s2)
        lhs = bool(t[0] and t[1])
        rhs = bool(t[0]) and bool(t[1])
        if lhs != rhs:
            all_cases_hold = False
            break
    print("This property is the definition of the logical AND operator, which holds for Python's `and` when converted to bool.")
    is_true = all_cases_hold
    print("Result: Tested all four truthiness combinations, and the equality holds. Statement F is True.")
    if is_true: true_statements.append('F')

    # --- Statement G ---
    print("\n--- Evaluating G ---")
    print("Statement: The expression (a and b) - (b and a) is always empty for any two sets a, b")
    a_g, b_g = {1}, {2}
    print(f"Counterexample: a = {a_g}, b = {b_g} (two non-empty, non-equal sets)")
    result = (a_g and b_g) - (b_g and a_g)
    print(f"(a and b) -> b. (b and a) -> a. Expression becomes b - a: {result}")
    is_true = (result == set())
    print(f"Result: The result is not an empty set. Statement G is False.")
    if is_true: true_statements.append('G')

    # --- Statement H ---
    print("\n--- Evaluating H ---")
    print("Statement: If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]")
    x_h = (False, 1)
    y_h = (True, 2)
    print(f"Counterexample: x = {x_h}, y = {y_h}")
    lhs = (x_h and y_h)[0]
    rhs = x_h[0] and y_h[0]
    print(f"x is truthy (non-empty tuple). 'x and y' -> y. LHS is y[0]: {lhs}")
    print(f"x[0] is False. 'x[0] and y[0]' -> x[0]. RHS is: {rhs}")
    is_true = (lhs == rhs)
    print(f"Result: {lhs} is not equal to {rhs}. Statement H is False.")
    if is_true: true_statements.append('H')

    # --- Statement I ---
    print("\n--- Evaluating I ---")
    print("Statement: For any sets p, q, r: (p and q) and r equals p and (q and r)")
    print("This checks the associativity of the 'and' operator. It evaluates to the first falsy operand, or the last if all are truthy.")
    truthy, falsy = [1], []
    is_true = True
    for p, q, r in itertools.product([truthy, falsy], repeat=3):
        if ((p and q) and r) is not (p and (q and r)):
            is_true = False
            break
    print("Result: Tested all 8 truthiness combinations, and the equality holds. Statement I is True.")
    if is_true: true_statements.append('I')

    # --- Statement J ---
    print("\n--- Evaluating J ---")
    print("Statement: all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b")
    a_j, b_j = {1, 2, 3}, {4, 5, 6}
    try:
        all(x and y for x, y in zip(a_j, b_j))
        print("Executing 'all(x and y for x, y in zip(a, b))' did not raise a TypeError.")
        is_true = False
    except TypeError:
        print("A TypeError was raised.")
        is_true = True
    print("Result: zip works on sets, and 'and' works on integers. No TypeError is raised. Statement J is False.")
    if is_true: true_statements.append('J')

    # --- Final Summary ---
    print("\n" + "="*40)
    print("Summary of True Statements:")
    true_statements.sort()
    print("The following statements are true:", ", ".join(true_statements))
    print("="*40)

if __name__ == '__main__':
    evaluate_statements()

<<<D, E, F, I>>>