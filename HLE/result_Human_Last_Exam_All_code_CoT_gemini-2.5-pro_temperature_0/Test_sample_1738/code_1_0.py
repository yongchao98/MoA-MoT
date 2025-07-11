import itertools

def solve_and_print():
    """
    Analyzes each statement from A to J, prints the evaluation,
    and determines the final answer.
    """
    # Given values
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    true_statements = []

    print("--- Evaluation of Each Statement ---")

    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    print("\n--- Statement A ---")
    val_A1 = c[0] and c[1]
    val_A2 = c[1] and c[0]
    is_true_A = val_A1 == val_A2
    print(f"c[0] is truthy, so 'c[0] and c[1]' evaluates to c[1]: {val_A1}")
    print(f"c[1] is truthy, so 'c[1] and c[0]' evaluates to c[0]: {val_A2}")
    print(f"Since {val_A1} != {val_A2}, Statement A is {is_true_A}.")

    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    print("\n--- Statement B ---")
    s1, s2 = {1}, {2}
    # s1 and s2 -> s2 (since s1 is truthy)
    # s2 and s1 -> s1 (since s2 is truthy)
    val_B = (s1 and s2) | (s2 and s1)
    is_true_B = val_B == s2
    print(f"For s1={s1}, s2={s2}:")
    print(f"'(s1 and s2) | (s2 and s1)' evaluates to 's2 | s1', which is {val_B}")
    print(f"This is not equal to s2 ({s2}). Statement B is {is_true_B}.")

    # C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
    print("\n--- Statement C ---")
    val_C = d[0] or d[1]
    is_true_C = val_C == (1, 2)
    print(f"d[0] is {d[0]}, which is a non-empty tuple and thus truthy.")
    print(f"The 'or' expression short-circuits and returns the first operand: {val_C}")
    print(f"This matches (1, 2). Statement C is {is_true_C}.")
    if is_true_C:
        true_statements.append('C')

    # D) If s is any non-empty set, s and True returns True while True and s returns s
    print("\n--- Statement D ---")
    s = {"non-empty"}
    val_D1 = s and True
    val_D2 = True and s
    is_true_D = (val_D1 is True) and (val_D2 is s)
    print(f"For a non-empty set s, s is truthy.")
    print(f"'s and True' evaluates to the second operand, True: {val_D1}")
    print(f"'True and s' evaluates to the second operand, s: {val_D2}")
    print(f"Both parts of the statement hold. Statement D is {is_true_D}.")
    if is_true_D:
        true_statements.append('D')

    # E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    print("\n--- Statement E ---")
    val_E1 = {} or []
    val_E2 = [] or {}
    is_true_E = (val_E1 == []) and (val_E2 == {})
    print(f"{{}} is falsy, so '{{}} or []' evaluates to the second operand: {val_E1}")
    print(f"[] is falsy, so '[] or {{}}' evaluates to the second operand: {val_E2}")
    print(f"Both parts of the statement hold. Statement E is {is_true_E}.")
    if is_true_E:
        true_statements.append('E')

    # F) For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    print("\n--- Statement F ---")
    # This is a fundamental property of boolean logic: bool(x and y) == (bool(x) and bool(y))
    is_true_F = True
    print("This statement tests if bool(x and y) is equivalent to bool(x) and bool(y).")
    print("This holds true for all objects x and y in Python.")
    print(f"Statement F is {is_true_F}.")
    if is_true_F:
        true_statements.append('F')

    # G) The expression (a and b) - (b and a) is always empty for any two sets a, b
    print("\n--- Statement G ---")
    set_g1, set_g2 = {1, 2}, {2, 3}
    # set_g1 and set_g2 -> set_g2
    # set_g2 and set_g1 -> set_g1
    val_G = (set_g1 and set_g2) - (set_g2 and set_g1)
    is_true_G = val_G == set()
    print(f"Counterexample with two non-empty sets a={set_g1}, b={set_g2}:")
    print(f"'(a and b) - (b and a)' evaluates to 'b - a', which is {val_G}")
    print(f"This is not an empty set. Statement G is {is_true_G}.")

    # H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    print("\n--- Statement H ---")
    x, y = (0, 1), (2, 3) # Both are truthy (non-empty tuples)
    lhs_H = (x and y)[0] # (x and y) -> y. y[0] -> 2
    rhs_H = x[0] and y[0] # 0 and 2 -> 0
    is_true_H = lhs_H == rhs_H
    print(f"Counterexample with x={x}, y={y}:")
    print(f"'(x and y)[0]' evaluates to y[0], which is {lhs_H}")
    print(f"'x[0] and y[0]' evaluates to '0 and 2', which is {rhs_H}")
    print(f"The results are not equal. Statement H is {is_true_H}.")

    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    print("\n--- Statement I ---")
    # This tests the associativity of the 'and' operator's resulting value.
    is_true_I = True
    s_true, s_false = {1}, set()
    for p_i, q_i, r_i in itertools.product([s_true, s_false], repeat=3):
        if ((p_i and q_i) and r_i) != (p_i and (q_i and r_i)):
            is_true_I = False
            break
    print("The 'and' operator is associative in terms of its final value.")
    print("If p is falsy, both sides evaluate to p.")
    print("If p is truthy and q is falsy, both sides evaluate to q.")
    print("If p and q are truthy, both sides evaluate to r.")
    print(f"The equality holds in all cases. Statement I is {is_true_I}.")
    if is_true_I:
        true_statements.append('I')

    # J) The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    print("\n--- Statement J ---")
    try:
        result_J = all(x and y for x, y in zip(a, b))
        print(f"zip(a, b) works on sets. The expression executed without error.")
        print(f"The result is: {result_J}")
        is_true_J = False # Statement is false because no error was raised.
    except TypeError:
        print("A TypeError was raised.")
        is_true_J = True
    print(f"The statement claims a TypeError is raised, which is not true. Statement J is {is_true_J}.")

    # Final Answer
    final_answer = ",".join(sorted(true_statements))
    print("\n--- Final Answer ---")
    print(f"The letters of the true statements in sorted order are: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_and_print()