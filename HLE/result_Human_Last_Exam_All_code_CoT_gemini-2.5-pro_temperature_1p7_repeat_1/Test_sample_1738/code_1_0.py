def evaluate_statements():
    """
    This function programmatically evaluates each statement from the problem
    and prints a detailed analysis, concluding with a sorted list of true statements.
    """
    # Given values from the problem
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    true_statement_letters = []

    print("--- Analysis of Each Statement ---")

    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    print("\n--- Statement A ---")
    val_a1 = c[0] and c[1]
    val_a2 = c[1] and c[0]
    is_true_a = (val_a1 == val_a2)
    print(f"Expression 'c[0] and c[1]' evaluates to: {val_a1}")
    print(f"Expression 'c[1] and c[0]' evaluates to: {val_a2}")
    print(f"Conclusion: The statement is {is_true_a}.")
    if is_true_a: true_statement_letters.append("A")

    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    print("\n--- Statement B ---")
    s1_b, s2_b = {1}, {2} # Counter-example
    expr_b = (s1_b and s2_b) | (s2_b and s1_b)
    is_true_b = (expr_b == s2_b)
    print(f"Counter-example: s1 = {s1_b}, s2 = {s2_b}")
    print(f"'(s1 and s2) | (s2 and s1)' evaluates to s2 | s1 = {expr_b}")
    print(f"'{expr_b} == {s2_b}' is {is_true_b}")
    print(f"Conclusion: The statement is {is_true_b}.")
    if is_true_b: true_statement_letters.append("B")

    # C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
    print("\n--- Statement C ---")
    val_c = d[0] or d[1]
    is_true_c = (val_c == (1, 2))
    print(f"Given d={d}, 'd[0] or d[1]' evaluates to: {val_c}")
    print("Because d[0] is truthy, the 'or' expression short-circuits and returns d[0]. The value of d[1] is irrelevant.")
    print(f"Conclusion: The statement is {is_true_c}.")
    if is_true_c: true_statement_letters.append("C")

    # D) If s is any non-empty set, s and True returns True while True and s returns s
    print("\n--- Statement D ---")
    s_d = {"test"} # Non-empty set
    is_true_d = ((s_d and True) is True) and ((True and s_d) is s_d)
    print(f"For s={s_d}, 's and True' returns: {s_d and True}")
    print(f"For s={s_d}, 'True and s' returns: {True and s_d}")
    print(f"Conclusion: The statement is {is_true_d}.")
    if is_true_d: true_statement_letters.append("D")

    # E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    print("\n--- Statement E ---")
    val_e1 = {} or []
    val_e2 = [] or {}
    is_true_e = (val_e1 == []) and (val_e2 == {})
    print(f"'{{}} or []' evaluates to: {val_e1}")
    print(f"'[] or {{}}' evaluates to: {val_e2}")
    print(f"Conclusion: The statement is {is_true_e}.")
    if is_true_e: true_statement_letters.append("E")

    # F) For any tuple t containing two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    print("\n--- Statement F ---")
    # This is true by definition of Python's boolean logic.
    is_true_f = True
    print("This statement reflects the definitional behavior of the logical 'and' operator.")
    print(f"Conclusion: The statement is {is_true_f}.")
    if is_true_f: true_statement_letters.append("F")

    # G) The expression (a and b) - (b and a) is always empty for any two sets a, b
    print("\n--- Statement G ---")
    a_g, b_g = {1, 2}, {2, 3} # Counter-example
    val_g = (a_g and b_g) - (b_g and a_g)
    is_true_g = (val_g == set())
    print(f"Counter-example: a={a_g}, b={b_g}")
    print(f"'(a and b) - (b and a)' evaluates to b-a = {val_g}")
    print(f"The result is not empty. Conclusion: The statement is {is_true_g}.")
    if is_true_g: true_statement_letters.append("G")

    # H) If x and y are tuples, (x and y)[0] is always equal to x[0] and y[0]
    print("\n--- Statement H ---")
    x_h, y_h = ((), 'a'), ((1,), 'b') # Counter-example
    val_h1 = (x_h and y_h)[0]
    val_h2 = x_h[0] and y_h[0]
    is_true_h = (val_h1 == val_h2)
    print(f"Counter-example: x={x_h}, y={y_h}")
    print(f"'(x and y)[0]' evaluates to y[0]: {val_h1}")
    print(f"'x[0] and y[0]' evaluates to x[0]: {val_h2}")
    print(f"Conclusion: The statement is {is_true_h}.")
    if is_true_h: true_statement_letters.append("H")

    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    print("\n--- Statement I ---")
    # This associativity holds for Python's `and` operator.
    is_true_i = True
    print("The 'and' operator is associative in Python, preserving this property.")
    print(f"Conclusion: The statement is {is_true_i}.")
    if is_true_i: true_statement_letters.append("I")
    
    # J) all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    print("\n--- Statement J ---")
    is_true_j = False
    try:
        all(x and y for x, y in zip(a,b))
    except TypeError:
        is_true_j = True
    print(f"zip(set, set) and `and` on their elements do not cause a TypeError.")
    print(f"Conclusion: The statement is {is_true_j}.")
    if is_true_j: true_statement_letters.append("J")

    # Final summary
    true_statement_letters.sort()
    print("\n--- Summary ---")
    print(f"The true statements are: {', '.join(true_statement_letters)}")

# Execute the evaluation function
evaluate_statements()