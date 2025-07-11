def verify_statements():
    """
    Provides a programmatic verification for each statement's truth value.
    Prints True/False corresponding to whether each statement is correct.
    """
    # --- Given values ---
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    print("--- Verification Results ---")

    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    # For truthy a, b: (a and b) -> b; (b and a) -> a. These are not equal.
    print(f"A: { (c[0] and c[1]) == (c[1] and c[0]) }")

    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    # Counterexample: s1={10}, s2={20}. Expression is (s2|s1) which is {10,20}, not s2.
    s1_b, s2_b = {10}, {20}
    print(f"B: { ((s1_b and s2_b) | (s2_b and s1_b)) == s2_b }")

    # C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
    # d[0] is truthy, so `or` short-circuits, returning d[0].
    print(f"C: { (d[0] or d[1]) == (1, 2) }")

    # D) If s is any non-empty set, s and True returns True while True and s returns s
    # A non-empty set s is truthy. (s and True) -> True. (True and s) -> s.
    s_d = a
    print(f"D: { (s_d and True is True) and (True and s_d is s_d) }")

    # E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    # Both {} and [] are falsy. `falsy or x` evaluates to `x`.
    print(f"E: { (({} or []) == []) and (([] or {}) == {}) }")

    # F) For any tuple t ... bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    # This is a logical identity.
    t_f = ({1}, set())
    print(f"F: { bool(t_f[0] and t_f[1]) == (bool(t_f[0]) and bool(t_f[1])) }")

    # G) The expression (a and b) - (b and a) is always empty for any two sets a, b
    # Counterexample: for non-empty a, b, this is `b - a`, which may not be empty.
    print(f"G: { ((a and b) - (b and a)) == set() }")

    # H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    # Counterexample: x=(0,1), y=(2,3). LHS is 2, RHS is 0.
    x_h, y_h = (0, 1), (2, 3)
    print(f"H: { (x_h and y_h)[0] == (x_h[0] and y_h[0]) }")

    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    # The `and` operator's evaluation is associative.
    p_i, q_i, r_i = {1}, set(), {3}
    print(f"I: { ((p_i and q_i) and r_i) == (p_i and (q_i and r_i)) }")

    # J) The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    # This code runs fine. So the statement that it raises an error is false.
    try:
        all(x and y for x, y in zip(a, b))
        # No error was raised, so the statement claiming an error occurs is False.
        print("J: False")
    except TypeError:
        # A TypeError occurred, so the statement is True.
        print("J: True")

verify_statements()