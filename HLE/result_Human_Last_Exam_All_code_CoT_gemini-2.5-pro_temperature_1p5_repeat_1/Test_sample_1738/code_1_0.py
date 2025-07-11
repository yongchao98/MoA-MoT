import collections.abc

def main():
    """
    Evaluates each statement from the problem description and prints the final answer.
    """
    # Given values
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    true_statements = []

    print("--- Evaluating Statements ---")

    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    print("\n--- A) The expression c[0] and c[1] always returns the same value as c[1] and c[0] ---")
    val1 = c[0] and c[1]
    val2 = c[1] and c[0]
    print(f"c[0] is {c[0]}, which is truthy.")
    print(f"c[1] is {c[1]}, which is truthy.")
    print(f"Therefore, 'c[0] and c[1]' evaluates to the second operand, c[1]: {val1}")
    print(f"And 'c[1] and c[0]' evaluates to the second operand, c[0]: {val2}")
    is_true_A = val1 == val2
    print(f"Statement A is {is_true_A}, because {val1} != {val2}.")
    if is_true_A: true_statements.append('A')

    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    print("\n--- B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2 ---")
    s1, s2 = {1}, {2}
    print(f"Let's test with a counter-example: s1 = {s1}, s2 = {s2}.")
    term1 = s1 and s2
    term2 = s2 and s1
    result = term1 | term2
    print(f"'s1 and s2' (s1 is truthy) evaluates to s2: {term1}")
    print(f"'s2 and s1' (s2 is truthy) evaluates to s1: {term2}")
    print(f"The full expression is (s2 | s1), which is the union: {result}")
    is_true_B = result == s2
    print(f"Statement B is {is_true_B}, because the union {result} is not equal to s2 {s2}.")
    if is_true_B: true_statements.append('B')
    
    # C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
    print("\n--- C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples ---")
    result_c = d[0] or d[1]
    print(f"For the given d, d[0] is {d[0]}, which is a non-empty tuple and thus truthy.")
    print("The 'or' operator short-circuits and returns the first truthy operand.")
    print(f"Therefore, 'd[0] or d[1]' evaluates to d[0]: {result_c}")
    is_true_C = result_c == (1, 2)
    print(f"Statement C is {is_true_C}.")
    if is_true_C: true_statements.append('C')

    # D) If s is any non-empty set, s and True returns True while True and s returns s
    print("\n--- D) If s is any non-empty set, s and True returns True while True and s returns s ---")
    s = {10, 20} # A non-empty set is truthy
    part1 = (s and True)
    part2 = (True and s)
    print(f"Let s = {s}. s is truthy.")
    print(f"'s and True' evaluates to the second operand, True: {part1}")
    print(f"'True and s' evaluates to the second operand, s: {part2}")
    is_true_D = (part1 is True) and (part2 is s)
    print(f"Statement D is {is_true_D}.")
    if is_true_D: true_statements.append('D')

    # E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    print("\n--- E) The expression {} or [] evaluates to [] while [] or {} evaluates to {} ---")
    part1 = {} or []
    part2 = [] or {}
    print(f"{} is falsey. '{{}} or []' evaluates to the second operand: {part1}")
    print(f"[] is falsey. '[] or {{}}' evaluates to the second operand: {part2}")
    is_true_E = (part1 == []) and (part2 == {})
    print(f"Statement E is {is_true_E}.")
    if is_true_E: true_statements.append('E')
    
    # F) For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    print("\n--- F) For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1]) ---")
    print("This is a fundamental property of how the 'and' operator's boolean result is determined.")
    t1, t2 = ({1}, {2}), ({}, {2})
    print(f"Test case 1 (both truthy): t = {t1}. LHS: bool({t1[0]} and {t1[1]}) -> bool({t1[1]}) -> True. RHS: bool({t1[0]}) and bool({t1[1]}) -> True and True -> True.")
    print(f"Test case 2 (first falsey): t = {t2}. LHS: bool({t2[0]} and {t2[1]}) -> bool({t2[0]}) -> False. RHS: bool({t2[0]}) and bool({t2[1]}) -> False and True -> False.")
    is_true_F = True
    print(f"Statement F is {is_true_F}.")
    if is_true_F: true_statements.append('F')
    
    # G) The expression (a and b) - (b and a) is always empty for any two sets a, b
    print("\n--- G) The expression (a and b) - (b and a) is always empty for any two sets a, b ---")
    print(f"Using the given non-empty sets a={a} and b={b}.")
    val1 = a and b
    val2 = b and a
    result = val1 - val2
    print(f"a is truthy, so 'a and b' evaluates to b: {val1}")
    print(f"b is truthy, so 'b and a' evaluates to a: {val2}")
    print(f"The expression becomes b - a, which is {result}.")
    is_true_G = result == set()
    print(f"Statement G is {is_true_G}, because the result is not an empty set.")
    if is_true_G: true_statements.append('G')

    # H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    print("\n--- H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0] ---")
    x, y = (0, 99), (8, 99)
    print(f"Let's test with a counter-example: x = {x}, y = {y}. Both are non-empty and truthy.")
    lhs = (x and y)[0]
    rhs = x[0] and y[0]
    print(f"LHS: '(x and y)[0]' -> y[0] -> {lhs}")
    print(f"RHS: 'x[0] and y[0]' -> '{x[0]} and {y[0]}' (0 is falsey) -> {rhs}")
    is_true_H = lhs == rhs
    print(f"Statement H is {is_true_H}, because {lhs} != {rhs}.")
    if is_true_H: true_statements.append('H')

    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    print("\n--- I) For any sets p, q, r: (p and q) and r equals p and (q and r) ---")
    print("This tests the associativity of the 'and' operator's return value. Let's trace the logic.")
    print("Case 1: p is falsey. LHS: (p and q)->p; (p) and r->p. RHS: p and (q and r)->p. They are equal.")
    print("Case 2: p is truthy. LHS: (p and q)->q; (q) and r->(q and r). RHS: p and (q and r)->(q and r). They are equal.")
    is_true_I = True
    print(f"The logic holds in all cases. Statement I is {is_true_I}.")
    if is_true_I: true_statements.append('I')

    # J) The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    print("\n--- J) The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b ---")
    try:
        list(zip(a,b)) # to show zip works
        print(f"zip(a, b) works on sets and produces pairs like: {list(zip(a,b))[:2]}...")
        all(x and y for x, y in zip(a, b))
        print("The code runs without raising a TypeError. The result of 'all' is True because all integer products are truthy.")
        is_true_J = False
    except TypeError:
        print("A TypeError was raised.")
        is_true_J = True
    print(f"Statement J is {is_true_J}.")
    if is_true_J: true_statements.append('J')
    
    # Final sorted answer
    final_answer = ", ".join(sorted(true_statements))
    print("\n--- Final Answer ---")
    print(f"The true statements are: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    main()