import sys

def py2_int_div(a, b):
    """Simulates Python 2's '/' operator for integers."""
    if isinstance(a, int) and isinstance(b, int):
        return a // b
    return a / b

def solve_and_print():
    """
    Analyzes and verifies statements about Python's division semantics.
    """
    true_statements = []

    print("Analyzing Python Division Semantics Statements:\n")

    # A) For any two positive integers x, y: x / y == x // y was true in Python 2 but false in Python 3
    print("--- Statement A ---")
    x, y = 5, 2
    py2_res = py2_int_div(x, y) == (x // y)
    py3_res = (x / y) == (x // y)
    print(f"Python 2 (simulated): {x} / {y} == {x} // {y} -> {py2_int_div(x, y)} == {x // y}, which is {py2_res}")
    print(f"Python 3: {x} / {y} == {x} // {y} -> {x / y} == {x // y}, which is {py3_res}")
    is_true = py2_res and not py3_res
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('A')

    # B) The expression 1/2 * 2 equals 1/2 + 1/2 in all Python versions
    print("--- Statement B ---")
    py2_lhs = py2_int_div(1, 2) * 2
    py2_rhs = py2_int_div(1, 2) + py2_int_div(1, 2)
    py3_lhs = 1/2 * 2
    py3_rhs = 1/2 + 1/2
    is_true = (py2_lhs == py2_rhs) and (py3_lhs == py3_rhs)
    print(f"Python 2 (simulated): ({py2_int_div(1, 2)}) * {2} == {py2_int_div(1, 2)} + {py2_int_div(1, 2)} -> {py2_lhs} == {py2_rhs}, which is {py2_lhs == py2_rhs}")
    print(f"Python 3: ({1}/{2}) * {2} == ({1}/{2}) + ({1}/{2}) -> {py3_lhs} == {py3_rhs}, which is {py3_lhs == py3_rhs}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('B')

    # C) For integers a, b: a / float(b) has always returned same result as float(a) / b
    print("--- Statement C ---")
    a, b = 7, 2
    res1 = a / float(b)
    res2 = float(a) / b
    is_true = res1 == res2
    print(f"In all versions, if one operand is a float, the other is promoted to float.")
    print(f"{a} / float({b}) == float({a}) / {b} -> {res1} == {res2}, which is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('C')

    # D) In Python 3, -7 // 3 and int(-7 / 3) always give same result
    print("--- Statement D ---")
    val1 = -7 // 3
    val2 = int(-7 / 3)
    is_true = val1 == val2
    print(f"-7 // 3 (floor division) evaluates to {val1}")
    print(f"int(-7 / 3) (truncation) evaluates to int({-7/3:.4f}) which is {val2}")
    print(f"Result: {val1} == {val2} is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('D')

    # E) For any non-zero integers x, y: x / y == float(x) / float(y) in Python 3
    print("--- Statement E ---")
    x, y = 5, 2
    res1 = x / y
    res2 = float(x) / float(y)
    is_true = res1 == res2
    print(f"In Python 3, '/' is true division, equivalent to float division.")
    print(f"{x} / {y} == float({x}) / float({y}) -> {res1} == {res2}, which is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('E')

    # F) The expression 5.0 // 2 has always returned a float across all versions
    print("--- Statement F ---")
    res = 5.0 // 2
    is_true = isinstance(res, float)
    print(f"The expression 5.0 // 2 evaluates to {res}.")
    print(f"The type of {res} is {type(res)}, so isinstance({res}, float) is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('F')

    # G) For positive integers n, m: (n//m)*m + n%m == n holds true in all versions
    print("--- Statement G ---")
    n, m = 17, 5
    res = (n // m) * m + n % m
    is_true = res == n
    print(f"This is the Division Algorithm Identity, which is fundamental in Python.")
    print(f"({n}//{m})*{m} + {n}%{m} == {n} -> ({n // m})*{m} + {n % m} == {n} -> {res} == {n}, which is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('G')

    # H) In Python 3, -3.0 // 2.0 == -(3.0 // 2.0)
    print("--- Statement H ---")
    lhs = -3.0 // 2.0
    rhs = -(3.0 // 2.0)
    is_true = lhs == rhs
    print(f"LHS: -3.0 // 2.0 -> floor of -1.5 is {lhs}")
    print(f"RHS: -(3.0 // 2.0) -> -(floor of 1.5) -> -({3.0 // 2.0}) is {rhs}")
    print(f"Result: {lhs} == {rhs} is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('H')

    # I) For any integers a, b: a / b in Python 3 is equivalent to a * (1.0/b)
    print("--- Statement I ---")
    a, b = 69, 10
    res1 = a / b
    res2 = a * (1.0 / b)
    is_true = res1 == res2
    print(f"This tests floating point precision differences.")
    print(f"{a} / {b} -> {res1}")
    print(f"{a} * (1.0/{b}) -> {a} * {1.0/b} -> {res2}")
    print(f"Result: {res1} == {res2} is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('I')

    # J) The expression sys.maxsize / 2 returns float in Python 3 even when result fits in int
    print("--- Statement J ---")
    res = sys.maxsize / 2
    is_true = isinstance(res, float)
    print(f"In Python 3, '/' always produces a float.")
    print(f"sys.maxsize / 2 results in {res}, which is of type {type(res)}.")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('J')

    # K) Floor division of negative numbers rounds toward negative infinity in all versions
    print("--- Statement K ---")
    val = -7 // 3
    is_true = (val == -3)
    print(f"This is definitional for Python's // operator.")
    print(f"-7 / 3 is approx -2.33. Rounding toward negative infinity gives -3.")
    print(f"Code check: -7 // 3 = {val}. The floor is indeed -3.")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('K')
    
    # L) For any integer n: (n/1) - n evaluates to 0 in Python 3 but their types may differ
    print("--- Statement L ---")
    n_large = 2**53 + 1
    res_large = (n_large/1.0) - n_large
    is_true = (res_large == 0)
    print(f"This fails for integers that lose precision when converted to float.")
    print(f"Example: n = 2**53 + 1 = {n_large}")
    print(f"(n/1.0) - n -> float({n_large}) - {n_large} -> {float(n_large)} - {n_large} = {res_large}")
    print(f"Result {res_large} != 0.")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('L')

    # M) Given integers a,b,c: if (a/b)/c == a/(b*c) is true numerically, it's also true for their types
    print("--- Statement M ---")
    a,b,c = 8,4,2
    type1 = type((a/b)/c)
    type2 = type(a/(b*c))
    is_true = type1 == type2
    print("In Python 3, both expressions always produce a float, so their types are always equal.")
    print(f"type(({a}/{b})/{c}) is {type1}")
    print(f"type({a}/({b}*{c})) is {type2}")
    print("Since the types are always equal, the implication is true.")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('M')

    # N) For positive integers x,y: sum([x/y]*y) equals x in Python 3 only if x is divisible by y
    print("--- Statement N ---")
    x1, y1 = 6, 3
    res1 = sum([x1/y1]*y1)
    if_part_holds = (res1 == x1)
    print(f"Test 'if' (divisible): sum([{x1}/{y1}]*{y1}) -> sum([{x1/y1}]*{y1}) -> {res1}. {res1} == {x1} is {if_part_holds}.")
    x2, y2 = 1, 3
    res2 = sum([x2/y2]*y2)
    only_if_part_holds = (res2 != x2)
    print(f"Test 'only if' (not divisible): sum([{x2}/{y2}]*{y2}) -> {res2:.17f}. {res2:.17f} == {x2} is {res2 == x2}.")
    is_true = if_part_holds and only_if_part_holds
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('N')

    # O) The expression (1/3) + (1/3) + (1/3) == 1.0 evaluates to True even with floating-point arithmetic
    print("--- Statement O ---")
    res = (1/3) + (1/3) + (1/3)
    is_true = res == 1.0
    print(f"This tests a classic floating point inaccuracy.")
    print(f"({1}/{3}) + ({1}/{3}) + ({1}/{3}) -> {1/3:.17f} + {1/3:.17f} + {1/3:.17f} = {res:.17f}")
    print(f"Result: {res:.17f} == 1.0 is {is_true}")
    print(f"Conclusion: The statement is {is_true}.\n")
    if is_true: true_statements.append('O')

    print("--- Summary ---")
    print("The following statements are true:")
    print(', '.join(true_statements))

solve_and_print()