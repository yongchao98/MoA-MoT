def solve_and_demonstrate():
    """
    Evaluates each statement from A to J, demonstrating with code whether
    it is true or false, and then prints the letters of the true statements.
    """
    # Initial given values
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))
    
    true_statements = []

    print("--- Analysis of Each Statement ---")

    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    val1_A = c[0] and c[1]
    val2_A = c[1] and c[0]
    print(f"\nA) Is `({c[0]}) and ({c[1]})` == `({c[1]}) and ({c[0]})`?")
    print(f"   LHS evaluates to: {val1_A}")
    print(f"   RHS evaluates to: {val2_A}")
    print(f"   Conclusion: False, because for truthy x and y, `x and y` is y, while `y and x` is x.")

    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    s1_B, s2_B = {1}, {2}
    result_B = (s1_B and s2_B) | (s2_B and s1_B)
    print(f"\nB) For s1={s1_B}, s2={s2_B}, is `(s1 and s2)|(s2 and s1)` == s2?")
    print(f"   `(s1 and s2)|(s2 and s1)` becomes `({s2_B})|({s1_B})` which evaluates to: {result_B}")
    print(f"   Conclusion: False, the expression results in the union s1|s2, not necessarily s2.")
    
    # C) The expression d[0] or d[1] evaluates to (1, 2)
    result_C = d[0] or d[1]
    print(f"\nC) Does `({d[0]}) or ({d[1]})` evaluate to (1, 2)?")
    print(f"   Expression evaluates to: {result_C}")
    print(f"   Conclusion: True. Since {d[0]} is truthy, the `or` expression short-circuits to the first operand.")
    true_statements.append('C')
    
    # D) If s is any non-empty set, s and True returns True while True and s returns s
    s_D = {'example'}
    result1_D = s_D and True
    result2_D = True and s_D
    print(f"\nD) For non-empty s={s_D}, is `s and True` == True, and `True and s` == s?")
    print(f"   `{s_D} and True` evaluates to: {result1_D}")
    print(f"   `True and {s_D}` evaluates to: {result2_D}")
    print(f"   Conclusion: True, this correctly describes `and` operator behavior.")
    true_statements.append('D')

    # E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    result1_E = {} or []
    result2_E = [] or {}
    print(f"\nE) Do `{{}} or []` == [] and `[] or {{}}` == {{}} hold?")
    print(f"   `{{}} or []` evaluates to: {result1_E}")
    print(f"   `[] or {{}}` evaluates to: {result2_E}")
    print(f"   Conclusion: True, as `falsey or x` evaluates to `x`.")
    true_statements.append('E')

    # F) For any tuple t containing two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    s1_F, s2_F = {1}, set() # A True/False example
    lhs_F = bool(s1_F and s2_F)
    rhs_F = bool(s1_F) and bool(s2_F)
    print(f"\nF) For s1={s1_F}, s2={s2_F}, is `bool(s1 and s2)` == `bool(s1) and bool(s2)`?")
    print(f"   LHS: `bool({s1_F} and {s2_F})` -> `bool({s1_F and s2_F})` -> {lhs_F}")
    print(f"   RHS: `bool({s1_F}) and bool({s2_F})` -> `{bool(s1_F)} and {bool(s2_F)}` -> {rhs_F}")
    print(f"   Conclusion: True, this is a fundamental property of boolean logic.")
    true_statements.append('F')
    
    # G) The expression (a and b) - (b and a) is always empty for any two sets a, b
    result_G = (a and b) - (b and a)
    print(f"\nG) For a={a}, b={b}, is `(a and b)-(b and a)` empty?")
    print(f"   Expression `({a} and {b}) - ({b} and {a})` becomes `{b} - {a}`, evaluating to: {result_G}")
    print(f"   Conclusion: False, this counterexample is not an empty set.")

    # H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    x_H, y_H = (0, 1), (2, 3)
    lhs_H = (x_H and y_H)[0]
    rhs_H = x_H[0] and y_H[0]
    print(f"\nH) For x={x_H}, y={y_H}, is `(x and y)[0]` == `x[0] and y[0]`?")
    print(f"   LHS `({x_H} and {y_H})[0]` -> `{y_H}[0]` -> {lhs_H}")
    print(f"   RHS `{x_H[0]} and {y_H[0]}` -> `{x_H[0]} and {y_H[0]}` -> {rhs_H}")
    print(f"   Conclusion: False, a counterexample where x[0] is falsey proves it.")
    
    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    p_I, q_I, r_I = {1}, set(), {3} # A True, False, True example
    lhs_I = (p_I and q_I) and r_I
    rhs_I = p_I and (q_I and r_I)
    print(f"\nI) For p={p_I}, q={q_I}, r={r_I}, is `(p and q) and r` == `p and (q and r)`?")
    print(f"   LHS `({p_I} and {q_I}) and {r_I}` -> `{p_I and q_I} and {r_I}` -> {lhs_I}")
    print(f"   RHS `{p_I} and ({q_I} and {r_I})` -> `{p_I} and {q_I and r_I}` -> {rhs_I}")
    print(f"   Conclusion: True, the `and` operator is associative.")
    true_statements.append('I')
    
    # J) The statement all(x and y for x, y in zip(a, b)) raises TypeError
    try:
        all(x and y for x, y in zip(a, b))
        result_J = "No TypeError raised."
    except TypeError as e:
        result_J = f"TypeError was raised: {e}"
    print(f"\nJ) Does `all(x and y for x,y in zip(a,b))` raise a TypeError?")
    print(f"   Evaluation result: {result_J}")
    print(f"   Conclusion: False, sets are iterable so `zip` works without error.")
    
    # Final Result
    print("\n" + "="*40)
    true_statements.sort()
    print("The sorted letter choices of all the true statements are:")
    print(",".join(true_statements))

solve_and_demonstrate()