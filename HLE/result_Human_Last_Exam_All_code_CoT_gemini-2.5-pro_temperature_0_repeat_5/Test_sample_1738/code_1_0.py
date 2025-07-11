def solve_semantics_question():
    """
    Evaluates each statement from A to J and prints the reasoning and result.
    Finally, it prints the sorted list of true statement letters.
    """
    # Given data
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    true_statements = []

    print("--- Evaluating Statement A ---")
    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    # `x and y` returns y if x is truthy. Both a and b are non-empty sets, so they are truthy.
    # `c[0] and c[1]` (i.e., `a and b`) evaluates to `b`.
    # `c[1] and c[0]` (i.e., `b and a`) evaluates to `a`.
    val1_A = c[0] and c[1]
    val2_A = c[1] and c[0]
    is_A_true = val1_A == val2_A
    print(f"The value of c[0] is: {c[0]}")
    print(f"The value of c[1] is: {c[1]}")
    print(f"The expression 'c[0] and c[1]' evaluates to: {val1_A}")
    print(f"The expression 'c[1] and c[0]' evaluates to: {val2_A}")
    print(f"Statement A is {is_A_true}.")
    if is_A_true: true_statements.append('A')
    print("-" * 30)

    print("--- Evaluating Statement B ---")
    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    # Let s1 = {1}, s2 = {2}. Both are non-empty and truthy.
    s1_B, s2_B = {1}, {2}
    # `s1_B and s2_B` evaluates to `s2_B`.
    # `s2_B and s1_B` evaluates to `s1_B`.
    # The expression becomes `s2_B | s1_B`, which is the union `{1, 2}`.
    # This is not equal to `s2_B`.
    result_B = (s1_B and s2_B) | (s2_B and s1_B)
    is_B_true = result_B == s2_B
    print(f"Using s1 = {s1_B}, s2 = {s2_B}")
    print(f"The expression '(s1 and s2) | (s2 and s1)' evaluates to: {result_B}")
    print(f"This is not equal to s2 ({s2_B}).")
    print(f"Statement B is {is_B_true}.")
    if is_B_true: true_statements.append('B')
    print("-" * 30)

    print("--- Evaluating Statement C ---")
    # C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
    # The claim "regardless of the values" makes this a universal statement.
    # Let's test a different tuple: d_alt = ((9, 9), (8, 8)).
    d_alt_C = ((9, 9), (8, 8))
    # `d_alt_C[0]` is truthy, so `d_alt_C[0] or d_alt_C[1]` evaluates to `d_alt_C[0]`.
    result_C = d_alt_C[0] or d_alt_C[1]
    # The result is (9, 9), not (1, 2). The counterexample proves the statement false.
    is_C_true = False
    print(f"For d_alt = {d_alt_C}, the expression 'd_alt[0] or d_alt[1]' evaluates to: {result_C}")
    print("This is not (1, 2), so the 'regardless' claim is false.")
    print(f"Statement C is {is_C_true}.")
    if is_C_true: true_statements.append('C')
    print("-" * 30)

    print("--- Evaluating Statement D ---")
    # D) If s is any non-empty set, s and True returns True while True and s returns s
    s_D = {1, 2, 3}  # A non-empty set is truthy.
    # `s_D and True`: `s_D` is truthy, so it returns the second operand, `True`.
    # `True and s_D`: `True` is truthy, so it returns the second operand, `s_D`.
    res1_D = s_D and True
    res2_D = True and s_D
    is_D_true = (res1_D is True) and (res2_D is s_D)
    print(f"Using s = {s_D}")
    print(f"The expression 's and True' evaluates to: {res1_D}")
    print(f"The expression 'True and s' evaluates to: {res2_D}")
    print(f"Statement D is {is_D_true}.")
    if is_D_true: true_statements.append('D')
    print("-" * 30)

    print("--- Evaluating Statement E ---")
    # E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    # `x or y` returns y if x is falsey. Both `{}` and `[]` are falsey.
    res1_E = {} or []
    res2_E = [] or {}
    is_E_true = (res1_E == []) and (res2_E == {})
    print(f"The expression '{{}} or []' evaluates to: {res1_E}")
    print(f"The expression '[] or {{}}' evaluates to: {res2_E}")
    print(f"Statement E is {is_E_true}.")
    if is_E_true: true_statements.append('E')
    print("-" * 30)

    print("--- Evaluating Statement F ---")
    # F) For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    # This is a fundamental property. The truthiness of `x and y` is always `bool(x) and bool(y)`.
    t = ({1}, {}) # Test case: truthy and falsey sets
    lhs_F = bool(t[0] and t[1]) # bool({1} and {}) -> bool({}) -> False
    rhs_F = bool(t[0]) and bool(t[1]) # bool({1}) and bool({}) -> True and False -> False
    is_F_true = lhs_F == rhs_F
    print(f"Using t = {t}")
    print(f"The expression 'bool(t[0] and t[1])' evaluates to: {lhs_F}")
    print(f"The expression 'bool(t[0]) and bool(t[1])' evaluates to: {rhs_F}")
    print(f"Statement F is {is_F_true}.")
    if is_F_true: true_statements.append('F')
    print("-" * 30)

    print("--- Evaluating Statement G ---")
    # G) The expression (a and b) - (b and a) is always empty for any two sets a, b
    # Using the given non-empty sets a and b.
    # `a and b` -> `b`. `b and a` -> `a`. The expression becomes `b - a`.
    result_G = (a and b) - (b and a)
    is_G_true = len(result_G) == 0
    print(f"Using a = {a}, b = {b}")
    print(f"The expression '(a and b) - (b and a)' evaluates to b - a, which is: {result_G}")
    print(f"The result is not empty.")
    print(f"Statement G is {is_G_true}.")
    if is_G_true: true_statements.append('G')
    print("-" * 30)

    print("--- Evaluating Statement H ---")
    # H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    # Counterexample: Let x[0] be falsey, e.g., 0.
    x_H, y_H = (0, 2), (3, 4)
    # `x_H` is a non-empty tuple, so it's truthy. `x_H and y_H` evaluates to `y_H`.
    lhs_H = (x_H and y_H)[0] # y_H[0] -> 3
    # `x_H[0]` is 0 (falsey). `x_H[0] and y_H[0]` evaluates to `x_H[0]`.
    rhs_H = x_H[0] and y_H[0] # 0 and 3 -> 0
    is_H_true = lhs_H == rhs_H
    print(f"Using x = {x_H}, y = {y_H}")
    print(f"The expression '(x and y)[0]' evaluates to: {lhs_H}")
    print(f"The expression 'x[0] and y[0]' evaluates to: {rhs_H}")
    print(f"Statement H is {is_H_true}.")
    if is_H_true: true_statements.append('H')
    print("-" * 30)

    print("--- Evaluating Statement I ---")
    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    # This tests the associativity of `and`. If p is falsey, both sides evaluate to p.
    # If p is truthy, both sides evaluate to `q and r`. The statement is true.
    p, q, r = {1}, {}, {3}
    lhs_I = (p and q) and r # ({1} and {}) -> {}; {} and {3} -> {}
    rhs_I = p and (q and r) # ({} and {3}) -> {}; {1} and {} -> {}
    is_I_true = lhs_I == rhs_I
    print(f"Using p={p}, q={q}, r={r}")
    print(f"The expression '(p and q) and r' evaluates to: {lhs_I}")
    print(f"The expression 'p and (q and r)' evaluates to: {rhs_I}")
    print(f"Statement I is {is_I_true}.")
    if is_I_true: true_statements.append('I')
    print("-" * 30)

    print("--- Evaluating Statement J ---")
    # J) The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    # `zip` works on sets (as iterables). `all` works on the resulting numbers. No TypeError occurs.
    is_J_true = False
    try:
        all(x and y for x, y in zip(a, b))
        print("The expression executed without a TypeError.")
    except TypeError:
        print("A TypeError was raised.")
        is_J_true = True
    print(f"Statement J is {is_J_true}.")
    if is_J_true: true_statements.append('J')
    print("-" * 30)

    # Final result
    final_answer = "".join(sorted(true_statements))
    print(f"\nThe true statements are: {', '.join(sorted(true_statements))}")
    print(f"\nFinal Answer in required format: <<<...>>>")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_semantics_question()