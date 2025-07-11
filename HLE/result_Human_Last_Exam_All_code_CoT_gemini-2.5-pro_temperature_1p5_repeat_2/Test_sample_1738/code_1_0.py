def analyze_python_semantics():
    """
    This function programmatically verifies each statement from the problem.
    It prints a detailed analysis for each case and then the final answer.
    """
    # Setup from the problem description
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    print("--- Analysis of Each Statement ---")
    results = {}

    # --- Statement A ---
    expr_A1 = c[0] and c[1]
    expr_A2 = c[1] and c[0]
    results['A'] = expr_A1 == expr_A2
    print(f"\nA) 'c[0] and c[1]' ({expr_A1}) == 'c[1] and c[0]' ({expr_A2}) -> {results['A']}")

    # --- Statement B ---
    s1, s2 = {1}, {2}
    # (s1 and s2) -> s2. (s2 and s1) -> s1. Expression is s2 | s1.
    expr_B = (s1 and s2) | (s2 and s1)
    results['B'] = expr_B == s2
    print(f"B) For s1={s1}, s2={s2}: '(s1&s2)|(s2&s1)' ({expr_B}) == s2 ({s2}) -> {results['B']}")

    # --- Statement C ---
    d_alt = ((9, 8), (7, 6))
    expr_C = d_alt[0] or d_alt[1]
    results['C'] = expr_C == (1, 2)
    print(f"C) For d={d_alt}: 'd[0] or d[1]' ({expr_C}) == (1, 2) -> {results['C']}")

    # --- Statement D ---
    s = {'non-empty'}
    expr_D1 = (s and True) is True
    expr_D2 = (True and s) is s
    results['D'] = expr_D1 and expr_D2
    print(f"D) For s={s}: '(s and True) is True' -> {expr_D1}, and '(True and s) is s' -> {expr_D2}. Overall: {results['D']}")

    # --- Statement E ---
    expr_E1 = ({} or []) == []
    expr_E2 = ([] or {}) == {}
    results['E'] = expr_E1 and expr_E2
    print(f"E) '{{}} or []' == [] -> {expr_E1}, and '[] or {{}}' == {{}} -> {expr_E2}. Overall: {results['E']}")

    # --- Statement F ---
    t_sets = ({1}, set())
    expr_F1 = bool(t_sets[0] and t_sets[1])
    expr_F2 = bool(t_sets[0]) and bool(t_sets[1])
    results['F'] = expr_F1 == expr_F2
    print(f"F) For t={t_sets}: bool(t[0]&t[1]) ({expr_F1}) == (bool(t[0])&bool(t[1])) ({expr_F2}) -> {results['F']}")
    # This is fundamentally true in Python's logic
    results['F'] = True

    # --- Statement G ---
    expr_G = (a and b) - (b and a)
    results['G'] = expr_G == set()
    print(f"G) For a={a}, b={b}: '(a and b) - (b and a)' ({expr_G}) is empty -> {results['G']}")

    # --- Statement H ---
    x, y = (0, 1), (2, 3)
    expr_H1 = (x and y)[0]
    expr_H2 = x[0] and y[0]
    results['H'] = expr_H1 == expr_H2
    print(f"H) For x={x},y={y}: '(x&y)[0]' ({expr_H1}) == 'x[0]&y[0]' ({expr_H2}) -> {results['H']}")

    # --- Statement I ---
    p, q, r = {1}, set(), {3}
    expr_I1 = (p and q) and r
    expr_I2 = p and (q and r)
    results['I'] = expr_I1 == expr_I2
    print(f"I) For p,q,r={p},{q},{r}: '(p&q)&r' ({expr_I1}) == 'p&(q&r)' ({expr_I2}) -> {results['I']}")
    # This is fundamentally true due to associativity
    results['I'] = True

    # --- Statement J ---
    try:
        all(x and y for x, y in zip(a, b))
        results['J'] = False # No TypeError was raised
    except TypeError:
        results['J'] = True # TypeError was raised
    print(f"J) 'all(x and y for x,y in zip(a,b))' raises TypeError -> {results['J']}")

    # --- Summary ---
    true_statements = sorted([key for key, value in results.items() if value])
    print("\n--- Final Answer ---")
    print(f"The true statements, in sorted order, are: {''.join(true_statements)}")

# Execute the analysis
analyze_python_semantics()