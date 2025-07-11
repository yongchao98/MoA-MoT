import sys

def analyze_division_statements():
    """
    Analyzes a series of statements about Python's division operators
    and prints a detailed evaluation for each.
    """
    results = {}
    print("Analyzing Python Division Semantics:")
    print("="*40)

    # A) For any two positive integers x, y: x / y == x // y was true in Python 2 but false in Python 3
    # Python 2 behavior: For positive integers, '/' was integer division, same as '//'.
    # E.g., in Python 2: 5 / 2 -> 2 and 5 // 2 -> 2. So the first part is true.
    # Python 3 behavior:
    x, y = 5, 2
    is_false_in_py3 = (x / y != x // y)
    explanation = f"True. In Py2, 5/2 is 2 and 5//2 is 2. In Py3, 5/2 is {x/y} and 5//2 is {x//y}, which are not equal."
    results['A'] = explanation if is_false_in_py3 else "False."

    # B) The expression 1/2 * 2 equals 1/2 + 1/2 in all Python versions
    # Python 2 behavior: (1/2)*2 -> 0*2 -> 0. (1/2)+(1/2) -> 0+0 -> 0. True in Py2.
    # Python 3 behavior:
    a = 1/2 * 2
    b = 1/2 + 1/2
    explanation = f"True. In Py3, (1/2 * 2) is {a} and (1/2 + 1/2) is {b}. In Py2, both expressions evaluate to 0."
    results['B'] = explanation if a == b else "False."

    # C) For integers a, b: a / float(b) has always returned same result as float(a) / b
    a, b = 10, 4
    res1 = a / float(b)
    res2 = float(a) / b
    explanation = f"True. Type promotion makes them equivalent. {a}/float({b}) is {res1} and float({a})/{b} is {res2}. The same holds for Py2."
    results['C'] = explanation if res1 == res2 else "False."

    # D) In Python 3, -7 // 3 and int(-7 / 3) always give same result
    a, b = -7, 3
    res1 = a // b
    res2 = int(a / b)
    explanation = f"False. -7 // 3 is {res1} (floor division). -7 / 3 is {a/b}, and int({a/b}) is {res2} (truncation)."
    results['D'] = explanation if res1 != res2 else "True."

    # E) For any non-zero integers x, y: x / y == float(x) / float(y) in Python 3
    x, y = 5, 2
    res1 = x / y
    res2 = float(x) / float(y)
    explanation = f"True. In Py3, '/' is always float division. {x}/{y} is {res1} and float({x})/float({y}) is {res2}."
    results['E'] = explanation if res1 == res2 else "False."

    # F) The expression 5.0 // 2 has always returned a float across all versions
    res = 5.0 // 2
    explanation = f"True. The result of 5.0 // 2 is {res} and its type is {type(res).__name__}. Floor division with a float operand returns a float."
    results['F'] = explanation if isinstance(res, float) else "False."

    # G) For positive integers n, m: (n//m)*m + n%m == n holds true in all versions
    n, m = 17, 5
    check = (n // m) * m + n % m
    explanation = f"True. This is the Division Algorithm Identity. For n={n}, m={m}: ({n//m})*{m} + {n%m} = {check}, which equals n."
    results['G'] = explanation if check == n else "False."

    # H) In Python 3, -3.0 // 2.0 == -(3.0 // 2.0)
    res1 = -3.0 // 2.0
    res2 = -(3.0 // 2.0)
    explanation = f"False. -3.0 // 2.0 evaluates to {res1}. -(3.0 // 2.0) evaluates to {res2}."
    results['H'] = explanation if res1 != res2 else "True."

    # I) For any integers a, b: a / b in Python 3 is equivalent to a * (1.0/b)
    a, b = 10, 3
    res1 = a / b
    res2 = a * (1.0 / b)
    explanation = f"True. Both expressions perform true division. {a}/{b} = {res1}. {a} * (1.0/{b}) = {res2}."
    results['I'] = explanation if res1 == res2 else "False."

    # J) The expression sys.maxsize / 2 returns float in Python 3 even when result fits in int
    res = sys.maxsize / 2
    explanation = f"True. The result of sys.maxsize / 2 is of type {type(res).__name__}, because '/' always returns a float in Py3."
    results['J'] = explanation if isinstance(res, float) else "False."

    # K) Floor division of negative numbers rounds toward negative infinity in all versions
    a, b = -5, 2
    res = a // b
    explanation = f"True. This is the definition of `//`. For example, {a} // {b} results in {res}."
    results['K'] = explanation if res == -3 else "False."

    # L) For any integer n: (n/1) - n evaluates to 0 in Python 3 but their types may differ
    n = 5
    res_val = (n / 1) - n
    eval_to_zero = (res_val == 0)
    types_differ = (type(res_val) is not type(0))
    explanation = f"True. For n={n}, ({n}/1)-{n} is {res_val} (type {type(res_val).__name__}), which equals 0. The type differs from int."
    results['L'] = explanation if eval_to_zero and types_differ else "False."

    # M) Given integers a,b,c: if (a/b)/c == a/(b*c) is true numerically, it's also true for their types
    a, b, c = 8, 2, 2
    res1 = (a / b) / c
    res2 = a / (b * c)
    types_equal = type(res1) == type(res2)
    explanation = f"True. In Py3, both expressions produce floats. type(({a}/{b})/{c}) is {type(res1).__name__}, type({a}/({b}*{c})) is {type(res2).__name__}."
    results['M'] = explanation if types_equal else "False."

    # N) For positive integers x,y: sum([x/y]*y) equals x in Python 3 only if x is divisible by y
    x1, y1 = 10, 2 # Divisible
    res1 = sum([x1 / y1] * y1)
    x2, y2 = 10, 3 # Not divisible
    res2 = sum([x2 / y2] * y2)
    cond1 = (res1 == x1)
    cond2 = (res2 != x2)
    explanation = f"True. Divisible case ({x1}, {y1}): sum = {res1} (=={x1}). Non-divisible case ({x2}, {y2}): sum = {res2} (!= {x2})."
    results['N'] = explanation if cond1 and cond2 else "False."

    # O) The expression (1/3) + (1/3) + (1/3) == 1.0 evaluates to True even with floating-point arithmetic
    res = (1/3) + (1/3) + (1/3)
    explanation = f"False. Due to floating point imprecision, (1/3)+(1/3)+(1/3) = {res}, which does not equal 1.0."
    results['O'] = explanation if res != 1.0 else "True."

    print("\n--- Evaluation Results ---\n")
    for statement, result_text in sorted(results.items()):
        print(f"[{statement}] {result_text}\n")
    
    true_statements = sorted([key for key, value in results.items() if not value.startswith('False')])
    print("="*40)
    print(f"Conclusion: The true statements are: {', '.join(true_statements)}")

analyze_division_statements()
