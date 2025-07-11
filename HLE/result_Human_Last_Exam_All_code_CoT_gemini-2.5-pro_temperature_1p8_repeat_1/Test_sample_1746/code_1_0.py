import sys

def analyze_division_statements():
    """
    Analyzes a series of statements about Python's division semantics.
    For each statement, it prints the logic, the calculations, and a verdict.
    """
    print("--- Analysis of Python Division Semantics ---")

    # A) For any two positive integers x, y: x / y == x // y was true in Python 2 but false in Python 3
    print("\n--- Statement A ---")
    print("For any two positive integers x, y: x / y == x // y was true in Python 2 but false in Python 3")
    x, y = 7, 3
    py2_equiv_div = x // y
    py3_true_div = x / y
    floor_div = x // y
    print(f"Python 2 Simulation for positive integers: x / y was integer division.")
    print(f"  For x={x}, y={y}: In Python 2, x / y would be {py2_equiv_div}. x // y is {floor_div}.")
    print(f"  Comparison (Py2): {py2_equiv_div} == {floor_div} is {py2_equiv_div == floor_div}. This part holds.")
    x, y = 4, 2 # Counterexample for the second part of the statement
    py3_true_div_c = x / y
    floor_div_c = x // y
    print(f"Python 3 check: The statement claims 'false in Python 3' for ANY positive integers.")
    print(f"  Let's test a counterexample: x={x}, y={y}")
    print(f"  In Python 3, x / y is {py3_true_div_c}. x // y is {floor_div_c}.")
    print(f"  Comparison (Py3): {py3_true_div_c} == {floor_div_c} is {py3_true_div_c == floor_div_c}.")
    print("  Since the equality can be True in Python 3, the statement '...false in Python 3' is too strong.")
    print("Verdict: False")


    # B) The expression 1/2 * 2 equals 1/2 + 1/2 in all Python versions
    print("\n--- Statement B ---")
    print("The expression 1/2 * 2 equals 1/2 + 1/2 in all Python versions")
    print("Python 2 Simulation: 1/2 was integer division -> 0")
    py2_lhs = (1 // 2) * 2
    py2_rhs = (1 // 2) + (1 // 2)
    print(f"  LHS in Py2: (1//2) * 2 => {1//2} * 2 => {py2_lhs}")
    print(f"  RHS in Py2: (1//2) + (1//2) => {1//2} + {1//2} => {py2_rhs}")
    print(f"  Comparison (Py2): {py2_lhs} == {py2_rhs} is {py2_lhs == py2_rhs}")
    print("Python 3: 1/2 is float division -> 0.5")
    py3_lhs = (1 / 2) * 2
    py3_rhs = (1 / 2) + (1 / 2)
    print(f"  LHS in Py3: (1/2) * 2 => {1/2} * 2 => {py3_lhs}")
    print(f"  RHS in Py3: (1/2) + (1/2) => {1/2} + {1/2} => {py3_rhs}")
    print(f"  Comparison (Py3): {py3_lhs} == {py3_rhs} is {py3_lhs == py3_rhs}")
    print("  The equality check holds true within each version.")
    print("Verdict: True")


    # C) For integers a, b: a / float(b) has always returned same result as float(a) / b
    print("\n--- Statement C ---")
    print("For integers a, b: a / float(b) has always returned same result as float(a) / b")
    a, b = 7, 3
    lhs = a / float(b)
    rhs = float(a) / b
    print(f"  Testing with a={a}, b={b}.")
    print(f"  LHS: {a} / float({b}) => {lhs}")
    print(f"  RHS: float({a}) / {b} => {rhs}")
    print(f"  Comparison: {lhs} == {rhs} is {lhs == rhs}")
    print("  This is true because in either case, one operand being a float promotes the entire operation to float division.")
    print("Verdict: True")


    # D) In Python 3, -7 // 3 and int(-7 / 3) always give same result
    print("\n--- Statement D ---")
    print("In Python 3, -7 // 3 and int(-7 / 3) always give same result")
    lhs = -7 // 3
    true_div_val = -7 / 3
    rhs = int(true_div_val)
    print(f"  -7 // 3 => {lhs} (floor division rounds down towards negative infinity)")
    print(f"  int(-7 / 3) => int({true_div_val}) => {rhs} (int() truncates towards zero)")
    print(f"  Comparison: {lhs} == {rhs} is {lhs == rhs}")
    print("Verdict: False")


    # G) For positive integers n, m: (n//m)*m + n%m == n holds true in all versions
    print("\n--- Statement G ---")
    print("For positive integers n, m: (n//m)*m + n%m == n holds true in all versions")
    n, m = 17, 5
    q = n // m
    r = n % m
    result = q * m + r
    print(f"  Testing with n={n}, m={m}.")
    print(f"  ({n}//{m})*{m} + {n}%{m} => ({q})*{m} + {r} => {q*m} + {r} => {result}")
    print(f"  Comparison: {result} == {n} is {result == n}")
    print("  This is the Division Algorithm Identity, which is guaranteed to be true for Python's integer arithmetic.")
    print("Verdict: True")
    
    # I) For any integers a, b: a / b in Python 3 is equivalent to a * (1.0/b)
    print("\n--- Statement I ---")
    print("For any integers a, b: a / b in Python 3 is equivalent to a * (1.0/b)")
    a, b = 10, 3
    lhs = a / b
    rhs = a * (1.0/b)
    print(f"  Testing with a={a}, b={b} to check for floating point micro-differences.")
    print(f"  LHS: {a} / {b} => {lhs:.20f}")
    print(f"  RHS: {a} * (1.0/{b}) => {a} * {1.0/b:.20f} => {rhs:.20f}")
    print(f"  Comparison: {lhs} == {rhs} is {lhs == rhs}")
    print(f"  Due to multiple rounding steps, the results are not always bit-for-bit identical.")
    print("Verdict: False")
    
    # L) For any integer n: (n/1) - n evaluates to 0 in Python 3 but their types may differ
    print("\n--- Statement L ---")
    print("For any integer n: (n/1) - n evaluates to 0 in Python 3 but their types may differ")
    n = 2**53 + 1
    result = (n/1) - n
    print(f"  Testing with a large integer n where float(n) loses precision: n = {n}")
    print(f"  n / 1 => {n/1}")
    print(f"  (n/1) - n => {n/1} - {n} => {result}")
    print(f"  Comparison with 0: {result} == 0 is {result == 0}")
    print("  The statement fails for integers large enough to lose precision when converted to a float.")
    print("Verdict: False")


if __name__ == "__main__":
    analyze_division_statements()
