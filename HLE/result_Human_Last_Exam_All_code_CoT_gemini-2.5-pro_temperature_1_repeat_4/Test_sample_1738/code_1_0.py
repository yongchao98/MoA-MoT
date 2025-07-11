def solve_and_print():
    """
    This function evaluates each statement and prints a detailed analysis.
    """
    # Given values
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))
    
    print("--- Analysis of Statements ---")

    # A)
    val_A1 = c[0] and c[1]
    val_A2 = c[1] and c[0]
    print(f"\nA) 'c[0] and c[1]' is {val_A1}, but 'c[1] and c[0]' is {val_A2}. Verdict: False")

    # B)
    s1, s2 = {1}, {2}
    val_B = (s1 and s2) | (s2 and s1)
    expected_B = s2
    print(f"\nB) For s1={s1}, s2={s2}, '(s1 and s2)|(s2 and s1)' is {val_B}, which is not {expected_B}. Verdict: False")

    # C)
    val_C = d[0] or d[1]
    expected_C = (1, 2)
    print(f"\nC) 'd[0] or d[1]' is {val_C}. Since d[0] is truthy, 'or' returns it. Verdict: True")
    
    # D)
    s = {100}
    val_D1 = s and True
    val_D2 = True and s
    print(f"\nD) For non-empty set s, 's and True' is {val_D1} and 'True and s' is {val_D2}. Verdict: True")

    # E)
    val_E1 = {} or []
    val_E2 = [] or {}
    print(f"\nE) '{{}} or []' is {repr(val_E1)} and '[] or {{}}' is {repr(val_E2)}. Verdict: True")

    # F)
    t_truthy = ({1}, {2})
    t_falsy = (set(), {2})
    val_F1 = bool(t_truthy[0] and t_truthy[1]) == (bool(t_truthy[0]) and bool(t_truthy[1]))
    val_F2 = bool(t_falsy[0] and t_falsy[1]) == (bool(t_falsy[0]) and bool(t_falsy[1]))
    print(f"\nF) The equality bool(x and y) == (bool(x) and bool(y)) holds for all x, y. Example results: {val_F1}, {val_F2}. Verdict: True")

    # G)
    val_G = (a and b) - (b and a)
    print(f"\nG) For non-empty a, b, '(a and b) - (b and a)' becomes 'b - a', which is {val_G}, not always empty. Verdict: False")

    # H)
    x, y = (0, 1), (2, 3)
    val_H1 = (x and y)[0]
    val_H2 = x[0] and y[0]
    print(f"\nH) For x={x}, y={y}, '(x and y)[0]' is {val_H1}, but 'x[0] and y[0]' is {val_H2}. Verdict: False")

    # I)
    p, q, r = {1}, {2}, {3}
    val_I1 = (p and q) and r
    val_I2 = p and (q and r)
    print(f"\nI) The 'and' operator is associative. E.g. for three truthy sets, ({p} and {q}) and {r} is {val_I1}, and {p} and ({q} and {r}) is {val_I2}. Verdict: True")
    
    # J)
    try:
        all(x and y for x, y in zip(a, b))
        result_J = "No TypeError"
    except TypeError:
        result_J = "TypeError raised"
    print(f"\nJ) 'zip' on sets works fine. The statement runs and produces '{result_J}'. Verdict: False")

    # Final result
    true_statements = ['C', 'D', 'E', 'F', 'I']
    print("\n--- Final Answer ---")
    print(f"The sorted list of true statements is: {','.join(sorted(true_statements))}")

solve_and_print()