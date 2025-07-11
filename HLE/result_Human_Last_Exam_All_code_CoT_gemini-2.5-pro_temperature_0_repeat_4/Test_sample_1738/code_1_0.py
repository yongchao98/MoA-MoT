def analyze_statements():
    """
    Analyzes each statement from A to J and prints the reasoning.
    Finally, it prints the sorted list of true statements.
    """
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    true_statements = []

    print("--- Analysis of Statements ---")

    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    print("\n--- Statement A ---")
    # c[0] is truthy, so `c[0] and c[1]` evaluates to c[1].
    # c[1] is truthy, so `c[1] and c[0]` evaluates to c[0].
    val_A1 = c[0] and c[1]
    val_A2 = c[1] and c[0]
    is_A_true = val_A1 == val_A2
    print(f"Given a = {a}, b = {b}, c = (a, b)")
    print(f"c[0] and c[1] -> a and b -> {val_A1}")
    print(f"c[1] and c[0] -> b and a -> {val_A2}")
    print(f"Result: {val_A1} == {val_A2} is {is_A_true}. Statement A is False.")

    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    print("\n--- Statement B ---")
    s1, s2 = {10}, {20} # Two non-empty sets
    # s1 is truthy, so `s1 and s2` -> s2
    # s2 is truthy, so `s2 and s1` -> s1
    # The expression becomes s2 | s1 (union)
    result_B = (s1 and s2) | (s2 and s1)
    is_B_true = result_B == s2
    print(f"Test with s1={s1}, s2={s2}")
    print(f"(s1 and s2) | (s2 and s1) -> ({s1 and s2}) | ({s2 and s1}) -> {result_B}")
    print(f"Result: {result_B} == {s2} is {is_B_true}. This is only true if s1 is a subset of s2. Statement B is False.")

    # C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
    print("\n--- Statement C ---")
    # d[0] is truthy, so `d[0] or d[1]` evaluates to d[0].
    # However, the "regardless of the values" part is key.
    d_counterexample = ((), (3, 4))
    result_C = d_counterexample[0] or d_counterexample[1]
    print(f"For the given d={d}, d[0] or d[1] -> {d[0] or d[1]}. This part is (1, 2).")
    print(f"But for a counterexample d={d_counterexample}, d[0] or d[1] -> {result_C}.")
    print(f"Result: The outcome depends on the values. Statement C is False.")

    # D) If s is any non-empty set, s and True returns True while True and s returns s
    print("\n--- Statement D ---")
    s = {"hello"} # A non-empty set is truthy
    res1_D = s and True
    res2_D = True and s
    is_D_true = (res1_D is True) and (res2_D is s)
    print(f"Test with s={s}")
    print(f"s and True -> {s} and True -> {res1_D}")
    print(f"True and s -> True and {s} -> {res2_D}")
    print(f"Result: The evaluations are correct. Statement D is True.")
    if is_D_true: true_statements.append("D")

    # E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    print("\n--- Statement E ---")
    # {} is falsy, so `{} or []` evaluates to [].
    # [] is falsy, so `[] or {}` evaluates to {}.
    res1_E = {} or []
    res2_E = [] or {}
    is_E_true = (res1_E == []) and (res2_E == {})
    print(f"{{}} or [] -> {res1_E}")
    print(f"[] or {{}} -> {res2_E}")
    print(f"Result: The evaluations are correct. Statement E is True.")
    if is_E_true: true_statements.append("E")

    # F) For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    print("\n--- Statement F ---")
    # This is a fundamental property of Python's boolean logic.
    print("bool(x and y) is defined to be logically equivalent to bool(x) and bool(y).")
    print("Result: This holds for any objects x and y, not just sets. Statement F is True.")
    true_statements.append("F")

    # G) The expression (a and b) - (b and a) is always empty for any two sets a, b
    print("\n--- Statement G ---")
    # For truthy a and b, this becomes b - a.
    result_G = (a and b) - (b and a)
    is_G_true = result_G == set()
    print(f"For a={a}, b={b}:")
    print(f"(a and b) - (b and a) -> ({a and b}) - ({b and a}) -> {result_G}")
    print(f"Result: The result is not an empty set. Statement G is False.")

    # H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    print("\n--- Statement H ---")
    x, y = (0, 1), (2, 3) # x[0] is falsy
    # x is truthy, so `x and y` -> y
    # `(x and y)[0]` -> y[0] -> 2
    # `x[0] and y[0]` -> 0 and 2 -> 0
    res1_H = (x and y)[0]
    res2_H = x[0] and y[0]
    is_H_true = res1_H == res2_H
    print(f"Test with x={x}, y={y}")
    print(f"(x and y)[0] -> ({x and y})[0] -> {res1_H}")
    print(f"x[0] and y[0] -> {x[0]} and {y[0]} -> {res2_H}")
    print(f"Result: {res1_H} == {res2_H} is {is_H_true}. Statement H is False.")

    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    print("\n--- Statement I ---")
    # This tests the associativity of the `and` operator's evaluation.
    p, q, r = {1}, {}, {3} # True, False, True
    res1_I = (p and q) and r
    res2_I = p and (q and r)
    is_I_true = res1_I == res2_I
    print(f"Test with p={p}, q={q}, r={r}")
    print(f"(p and q) and r -> (({p} and {q})) and {r} -> ({p and q}) and {r} -> {res1_I}")
    print(f"p and (q and r) -> {p} and (({q} and {r})) -> {p} and ({q and r}) -> {res2_I}")
    print(f"Result: The `and` operator is associative in its evaluation. Statement I is True.")
    if is_I_true: true_statements.append("I")

    # J) The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    print("\n--- Statement J ---")
    try:
        result_J = all(x and y for x, y in zip(a, b))
        print(f"Executing `all(x and y for x, y in zip(a, b))` succeeded and returned {result_J}.")
        print("Result: No TypeError was raised. Statement J is False.")
    except TypeError:
        print("Result: A TypeError was raised. Statement J is True.")

    print("\n--- Conclusion ---")
    final_answer = ", ".join(sorted(true_statements))
    print(f"The true statements are: {final_answer}")
    return final_answer

if __name__ == '__main__':
    analyze_statements()