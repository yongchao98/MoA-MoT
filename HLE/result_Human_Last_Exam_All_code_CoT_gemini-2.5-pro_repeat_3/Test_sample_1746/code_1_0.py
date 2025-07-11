import sys
import math

def analyze_statements():
    """
    Analyzes statements about Python's division semantics.
    """
    true_statements = []

    print("--- Analyzing Python Division Statements ---")

    # --- Statement A ---
    print("\n--- Statement A ---")
    print("A) For any two positive integers x, y: x / y == x // y was true in Python 2 but false in Python 3")
    # In Python 2, for integers, '/' was floor division, same as '//'. So `5/2` was 2 and `5//2` was 2. The first part is true.
    # Let's test the second part in Python 3.
    x1, y1 = 5, 2
    res_a1 = (x1 / y1 == x1 // y1)
    print(f"In Python 3, for x={x1}, y={y1}: ({x1}/{y1}) == ({x1}//{y1}) is ({x1/y1}) == ({x1//y1}), which is {res_a1}.")
    x2, y2 = 4, 2
    res_a2 = (x2 / y2 == x2 // y2)
    print(f"In Python 3, for x={x2}, y={y2}: ({x2}/{y2}) == ({x2}//{y2}) is ({x2/y2}) == ({x2//y2}), which is {res_a2}.")
    print("Conclusion: The statement claims it's always false in Python 3, but for x=4, y=2, it's True. Therefore, statement A is False.")

    # --- Statement B ---
    print("\n--- Statement B ---")
    print("B) The expression 1/2 * 2 equals 1/2 + 1/2 in all Python versions")
    print("In Python 2 (integer division): (1/2)*2 -> 0*2 -> 0. And (1/2)+(1/2) -> 0+0 -> 0. So 0 == 0 is True.")
    lhs_b = 1/2 * 2
    rhs_b = 1/2 + 1/2
    print(f"In Python 3 (float division): (1/2)*2 is {lhs_b}. And (1/2)+(1/2) is {rhs_b}. So {lhs_b} == {rhs_b} is {lhs_b == rhs_b}.")
    print("Conclusion: The equality holds in both Python 2 and Python 3. Therefore, statement B is True.")
    true_statements.append("B")

    # --- Statement C ---
    print("\n--- Statement C ---")
    print("C) For integers a, b: a / float(b) has always returned same result as float(a) / b")
    a, b = 7, 2
    res1_c = a / float(b)
    res2_c = float(a) / b
    print(f"For a={a}, b={b}: a / float(b) is {a} / {float(b)} = {res1_c}. float(a) / b is {float(a)} / {b} = {res2_c}.")
    print("In both Python 2 and 3, when one operand is a float, the other is promoted to a float and float division is performed. The results are identical.")
    print("Conclusion: Statement C is True.")
    true_statements.append("C")

    # --- Statement D ---
    print("\n--- Statement D ---")
    print("D) In Python 3, -7 // 3 and int(-7 / 3) always give same result")
    val1_d = -7 // 3
    val2_d = int(-7 / 3)
    print(f"In Python 3: -7 // 3 is {val1_d}. (-7/3 is {-7/3}).")
    print(f"int(-7 / 3) is int({-7/3}) which is {val2_d}.")
    print(f"The results {val1_d} and {val2_d} are not the same.")
    print("Conclusion: `//` rounds to negative infinity, `int()` truncates toward zero. Therefore, statement D is False.")

    # --- Statement E ---
    print("\n--- Statement E ---")
    print("E) For any non-zero integers x, y: x / y == float(x) / float(y) in Python 3")
    x, y = 10, 4
    res1_e = x / y
    res2_e = float(x) / float(y)
    print(f"In Python 3, for x={x}, y={y}: x / y is {x/y}. float(x) / float(y) is {float(x)/float(y)}. They are equal: {res1_e == res2_e}.")
    print("This is because the `/` operator in Python 3 performs 'true' division, producing a float, which is the same as dividing two floats.")
    print("Conclusion: Statement E is True.")
    true_statements.append("E")

    # --- Statement F ---
    print("\n--- Statement F ---")
    print("F) The expression 5.0 // 2 has always returned a float across all versions")
    res_f = 5.0 // 2
    print(f"In Python 3: 5.0 // 2 is {res_f}, and its type is {type(res_f)}.")
    print("In Python 2, 5.0 // 2 also resulted in 2.0, a float. When an operand in `//` is a float, the result is a float.")
    print("Conclusion: Statement F is True.")
    true_statements.append("F")

    # --- Statement G ---
    print("\n--- Statement G ---")
    print("G) For positive integers n, m: (n//m)*m + n%m == n holds true in all versions")
    n, m = 17, 5
    res_g = (n//m)*m + n%m
    print(f"For n={n}, m={m}: (n//m)*m + n%m is ({n//m})*{m} + {n%m} = {res_g}. The result {res_g} equals n ({n}).")
    print("This is the definition of the Division Identity, which Python's `//` and `%` operators are designed to satisfy for integers.")
    print("Conclusion: Statement G is True.")
    true_statements.append("G")
    
    # --- Statement H ---
    print("\n--- Statement H ---")
    print("H) In Python 3, -3.0 // 2.0 == -(3.0 // 2.0)")
    lhs_h = -3.0 // 2.0
    rhs_h = -(3.0 // 2.0)
    print(f"LHS: -3.0 // 2.0 (floor of -1.5) is {lhs_h}.")
    print(f"RHS: -(3.0 // 2.0) (-(floor of 1.5)) is -(1.0) which is {rhs_h}.")
    print(f"The results {lhs_h} and {rhs_h} are not equal.")
    print("Conclusion: Statement H is False.")

    # --- Statement I ---
    print("\n--- Statement I ---")
    print("I) For any integers a, b: a / b in Python 3 is equivalent to a * (1.0/b)")
    a, b = 10, 4
    res1_i = a / b
    res2_i = a * (1.0/b)
    print(f"For a={a}, b={b}: a / b is {res1_i}. a * (1.0/b) is {a} * ({1.0/b}) = {res2_i}. Equal: {res1_i == res2_i}.")
    print("In Python 3, `a/b` is float division, which is mathematically equivalent to `a * (1/b)`.")
    print("Conclusion: Statement I is True.")
    true_statements.append("I")

    # --- Statement J ---
    print("\n--- Statement J ---")
    print("J) The expression sys.maxsize / 2 returns float in Python 3 even when result fits in int")
    res_j = sys.maxsize / 2
    print(f"sys.maxsize / 2 results in a value of type {type(res_j)}.")
    print("In Python 3, the `/` operator always returns a float, regardless of the inputs or if the result is a whole number.")
    print("Conclusion: Statement J is True.")
    true_statements.append("J")

    # --- Statement K ---
    print("\n--- Statement K ---")
    print("K) Floor division of negative numbers rounds toward negative infinity in all versions")
    val_k = -7 // 3
    print(f"The true value of -7 / 3 is {-7/3}. Rounding toward negative infinity gives {math.floor(-7/3)}.")
    print(f"The result of -7 // 3 in Python is {val_k}. This behavior is consistent across versions.")
    print("Conclusion: Statement K is True.")
    true_statements.append("K")

    # --- Statement L ---
    print("\n--- Statement L ---")
    print("L) For any integer n: (n/1) - n evaluates to 0 in Python 3 but their types may differ")
    n = 5
    res_l = (n/1) - n
    print(f"For n={n}: (n/1) - n is ({n/1}) - {n} = {res_l}. Type is {type(res_l)}.")
    print(f"Does the result evaluate to 0? {res_l == 0}. Is its type different from type(0)? {type(res_l) != type(0)}.")
    print("Conclusion: Statement L is True.")
    true_statements.append("L")
    
    # --- Statement M ---
    print("\n--- Statement M ---")
    print("M) Given integers a,b,c: if (a/b)/c == a/(b*c) is true numerically, it's also true for their types")
    a,b,c = 10,2,2
    lhs_m = (a/b)/c
    rhs_m = a/(b*c)
    print(f"In Python 3, for a,b,c={a},{b},{c}:")
    print(f"LHS type: type((a/b)/c) is {type(lhs_m)}. RHS type: type(a/(b*c)) is {type(rhs_m)}.")
    print("In Python 3, both expressions always produce a float. In Python 2, both would produce an int. Thus, if they are numerically equal, their types are also always equal.")
    print("Conclusion: Statement M is True.")
    true_statements.append("M")

    # --- Statement N ---
    print("\n--- Statement N ---")
    print("N) For positive integers x,y: sum([x/y]*y) equals x in Python 3 only if x is divisible by y")
    x, y = 7, 3
    res_n = sum([x/y]*y)
    print(f"Test case where x is NOT divisible by y: x={x}, y={y}.")
    print(f"sum([x/y]*y) is sum([{x/y}]*{y}) = {res_n}. Does this equal x? {res_n == x}.")
    print(f"The equality {res_n} == {x} holds, but {x} is not divisible by {y}. This contradicts the 'only if' condition.")
    print("Conclusion: Statement N is False.")

    # --- Statement O ---
    print("\n--- Statement O ---")
    print("O) The expression (1/3) + (1/3) + (1/3) == 1.0 evaluates to True even with floating-point arithmetic")
    print("In Python 2, this is 0 + 0 + 0 == 1.0, which is False.")
    res_o = (1/3) + (1/3) + (1/3)
    print(f"In Python 3, (1/3) + (1/3) + (1/3) is {res_o}. The expression {res_o} == 1.0 is {res_o == 1.0}.")
    print("Since the statement is not qualified by a Python version and is False in Python 2, it is not universally true.")
    print("Conclusion: Statement O is False.")

    print("\n--- Summary ---")
    print(f"The true statements are: {', '.join(true_statements)}")


if __name__ == '__main__':
    analyze_statements()
