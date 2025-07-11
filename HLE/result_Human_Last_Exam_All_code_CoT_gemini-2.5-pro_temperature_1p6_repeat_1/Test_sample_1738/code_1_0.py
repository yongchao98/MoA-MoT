def evaluate_statements():
    """
    Analyzes Python semantic statements and prints the letters of the true ones.
    """
    # Given values for context in some questions
    a_given = {1, 2, 3}
    b_given = {3, 4, 5}
    c_given = (a_given, b_given)
    d_given = ((1, 2), (3, 4))

    true_statements = []

    # --- Statement A ---
    # The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    # For truthy a and b, (a and b) -> b, while (b and a) -> a. Since a != b, this is false.
    if (c_given[0] and c_given[1]) != (c_given[1] and c_given[0]):
        pass  # Statement A is False
    else:
        true_statements.append('A')

    # --- Statement B ---
    # For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    # For non-empty sets, this simplifies to s2 | s1, the union, which doesn't equal s2 in general.
    s1_b, s2_b = {1}, {2}
    if ((s1_b and s2_b) | (s2_b and s1_b)) != s2_b:
        pass  # Statement B is False
    else:
        true_statements.append('B')

    # --- Statement C ---
    # The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
    # The "regardless" part is key. Let's test a case where d[0] is falsy.
    d_test = ((), (3, 4))
    if (d_test[0] or d_test[1]) != (1, 2):
        pass # Statement C is False because the outcome depends on the values.
    else:
        true_statements.append('C')
        
    # --- Statement D ---
    # If s is any non-empty set, s and True returns True while True and s returns s
    # A non-empty set is truthy. `truthy and True` -> `True`. `True and truthy` -> `truthy`. This holds.
    s_d = {"test"}
    if (s_d and True is True) and (True and s_d is s_d):
        true_statements.append('D')

    # --- Statement E ---
    # The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    # `{}` is falsy, so `{} or []` -> `[]`. `[]` is falsy, so `[] or {}` -> `{}`. This is correct.
    if (({} or []) == []) and (([] or {}) == {}):
        true_statements.append('E')
        
    # --- Statement F ---
    # For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    # This logical equivalence holds for any two values, not just sets.
    # Case 1 (T, T): bool({1} and {2}) == bool({1}) and bool({2}) -> bool({2}) == True and True -> True == True
    # Case 2 (T, F): bool({1} and {}) == bool({1}) and bool({}) -> bool({}) == True and False -> False == False
    # Case 3 (F, T): bool({} and {1}) == bool({}) and bool({1}) -> bool({}) == False and True -> False == False
    # Case 4 (F, F): bool({} and {}) == bool({}) and bool({}) -> bool({}) == False and False -> False == False
    true_statements.append('F')

    # --- Statement G ---
    # The expression (a and b) - (b and a) is always empty for any two sets a, b
    # For non-empty a, b, this is `b - a`, which is not always empty.
    a_g, b_g = {1, 2}, {2, 3}
    if ((a_g and b_g) - (b_g and a_g)) != set():
        pass # Statement G is False
    else:
        true_statements.append('G')

    # --- Statement H ---
    # If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    # Counterexample: x=(False, 1), y=(True, 2).
    # (x and y)[0] -> y[0] -> True.
    # x[0] and y[0] -> False and True -> False. Not equal.
    x_h, y_h = (False, 1), (True, 2)
    if (x_h and y_h)[0] != (x_h[0] and y_h[0]):
        pass # Statement H is False
    else:
        true_statements.append('H')
        
    # --- Statement I ---
    # For any sets p, q, r: (p and q) and r equals p and (q and r)
    # The `and` operation is associative in this manner. It finds the first falsy operand or the last operand.
    # (p falsy): (p and q)->p; (p and r)->p.  ||  p and (...)->p.  Result: p == p.
    # (p truthy, q falsy): (p and q)->q; (q and r)->q.  ||  (q and r)->q; p and q -> q. Result: q == q.
    # ...and so on. The logic holds.
    true_statements.append('I')

    # --- Statement J ---
    # The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    # Sets are iterable, so zip works fine. The `and` operator works on any objects. `all` takes an iterable. No TypeError occurs.
    try:
        all(x and y for x, y in zip(a_given, b_given))
        pass # No TypeError was raised, so the statement is False
    except TypeError:
        true_statements.append('J') # This branch is not taken

    print("".join(sorted(true_statements)))

evaluate_statements()