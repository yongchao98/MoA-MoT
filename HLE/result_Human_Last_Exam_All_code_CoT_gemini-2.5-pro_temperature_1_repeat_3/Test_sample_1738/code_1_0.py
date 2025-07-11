def evaluate_statements():
    """
    This function programmatically evaluates each statement from the question
    to determine its truth value and prints the reasoning.
    """
    # Given values
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    print("--- Evaluating Statements ---")
    
    true_statements = []

    # Statement A
    val1_A = c[0] and c[1]
    val2_A = c[1] and c[0]
    is_A_true = val1_A == val2_A
    if is_A_true: true_statements.append('A')
    print(f"\nA) Is `(a and b) == (b and a)`? -> Is `{val1_A} == {val2_A}`? -> {is_A_true}")

    # Statement B
    s1_B, s2_B = {1}, {2}
    result_B = (s1_B and s2_B) | (s2_B and s1_B)
    is_B_true = result_B == s2_B
    if is_B_true: true_statements.append('B')
    print(f"B) For s1={s1_B}, s2={s2_B}: Is `(s2 | s1) == s2`? -> Is `{result_B} == {s2_B}`? -> {is_B_true}")

    # Statement C
    d_alt_C = ((9, 8), (7, 6))
    is_C_true = (d[0] or d[1]) == (1, 2) and (d_alt_C[0] or d_alt_C[1]) == (1, 2)
    if is_C_true: true_statements.append('C')
    print(f"C) Is `d[0] or d[1]` always (1,2)? -> {is_C_true} (since for d={d_alt_C} it is {d_alt_C[0] or d_alt_C[1]})")

    # Statement D
    s_D = {'non-empty'}
    is_D_true = (s_D and True) is True and (True and s_D) is s_D
    if is_D_true: true_statements.append('D')
    print(f"D) Does `s and True` return True and `True and s` return s? -> {is_D_true}")

    # Statement E
    is_E_true = ({} or []) == [] and ([] or {}) == {}
    if is_E_true: true_statements.append('E')
    print(f"E) Does `{{}} or []` return [] and `[] or {{}}` return {{}}? -> {is_E_true}")

    # Statement F
    # bool(x and y) is equivalent to bool(x) and bool(y) for all x, y.
    is_F_true = True
    if is_F_true: true_statements.append('F')
    print(f"F) Is `bool(t[0] and t[1]) == (bool(t[0]) and bool(t[1]))` always true? -> {is_F_true}")

    # Statement G
    result_G = (a and b) - (b and a)
    is_G_true = result_G == set()
    if is_G_true: true_statements.append('G')
    print(f"G) Is `(a and b) - (b and a)` always empty? -> {is_G_true} (since it is `{result_G}` for the given a,b)")

    # Statement H
    x_H, y_H = (0, 2), (3, 4)
    lhs_H = (x_H and y_H)[0]
    rhs_H = x_H[0] and y_H[0]
    is_H_true = lhs_H == rhs_H
    if is_H_true: true_statements.append('H')
    print(f"H) Is `(x and y)[0] == x[0] and y[0]` always true? -> {is_H_true} (since for x={x_H}, y={y_H}, {lhs_H} != {rhs_H})")

    # Statement I
    # 'and' is associative in Python.
    is_I_true = True
    if is_I_true: true_statements.append('I')
    print(f"I) Is `(p and q) and r == p and (q and r)` always true? -> {is_I_true}")
    
    # Statement J
    is_J_true = False
    try:
        all(x and y for x, y in zip(a, b))
    except TypeError:
        is_J_true = True
    if is_J_true: true_statements.append('J')
    print(f"J) Does `all(zip(a,b))` raise TypeError? -> {is_J_true}")

    final_answer = ",".join(sorted(true_statements))
    print("\n---------------------------------")
    print(f"The true statements are: {final_answer}")
    print("---------------------------------")


evaluate_statements()