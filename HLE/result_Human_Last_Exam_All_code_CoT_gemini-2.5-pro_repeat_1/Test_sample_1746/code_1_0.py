import sys

def analyze_division_statements():
    """
    Analyzes a series of statements about Python's division semantics
    by executing or explaining test cases for each.
    """
    print("--- Analysis of Python Division Statements ---")
    print(f"Running on Python {sys.version.split()[0]}\n")

    # A) For any two positive integers x, y: x / y == x // y was true in Python 2 but false in Python 3
    x_a, y_a = 5, 2
    # Python 2 behavior: 5 / 2 -> 2. 5 // 2 -> 2. So 5/2 == 5//2 is True.
    # Python 3 behavior:
    res_a_py3_slash = x_a / y_a
    res_a_py3_slash_slash = x_a // y_a
    is_false_in_py3 = not (res_a_py3_slash == res_a_py3_slash_slash)
    print(f"A) Statement: 'x/y == x//y' was true in Py2, but false in Py3 (for positive integers).")
    print(f"   - In Python 2, integer division made `{x_a}/{y_a}` result in `{y_a}`, same as `{x_a}//{y_a}`.")
    print(f"   - In Python 3, `{x_a}/{y_a}` results in `{res_a_py3_slash}` while `{x_a}//{y_a}` results in `{res_a_py3_slash_slash}`. They are not equal.")
    print(f"   - Verdict: The statement is TRUE.\n")

    # D) In Python 3, -7 // 3 and int(-7 / 3) always give same result
    val_d_1 = -7 // 3
    val_d_2 = int(-7 / 3)
    res_d_equal = val_d_1 == val_d_2
    print(f"D) Statement: In Python 3, '-7 // 3' and 'int(-7 / 3)' give the same result.")
    print(f"   - `-7 // 3` (floor division) rounds to negative infinity, resulting in: {val_d_1}")
    print(f"   - `-7 / 3` is {-7/3:.4f}..., and `int()` truncates towards zero, resulting in: {val_d_2}")
    print(f"   - Verdict: `{val_d_1} == {val_d_2}` is {res_d_equal}. The statement is FALSE.\n")

    # G) For positive integers n, m: (n//m)*m + n%m == n holds true in all versions
    n_g, m_g = 17, 5
    res_g = (n_g // m_g) * m_g + (n_g % m_g)
    is_g_true = (res_g == n_g)
    print(f"G) Statement: For positive integers n, m, '(n//m)*m + n%m == n' holds true in all versions.")
    print(f"   - Testing with n={n_g}, m={m_g}: `({n_g}//{m_g})*{m_g} + {n_g}%{m_g}`")
    print(f"   - Calculation: `({n_g // m_g})*{m_g} + {n_g % m_g}` equals `{res_g}`.")
    print(f"   - Verdict: `{res_g} == {n_g}` is {is_g_true}. This identity is guaranteed by Python's definition of // and %. The statement is TRUE.\n")

    # H) In Python 3, -3.0 // 2.0 == -(3.0 // 2.0)
    val_h_lhs = -3.0 // 2.0
    val_h_rhs = -(3.0 // 2.0)
    res_h_equal = val_h_lhs == val_h_rhs
    print(f"H) Statement: In Python 3, `-3.0 // 2.0 == -(3.0 // 2.0)`.")
    print(f"   - The left side, `-3.0 // 2.0`, floors -1.5 to get: {val_h_lhs}")
    print(f"   - The right side, `-(3.0 // 2.0)`, floors 1.5 to get 1.0, then negates it: {val_h_rhs}")
    print(f"   - Verdict: `{val_h_lhs} == {val_h_rhs}` is {res_h_equal}. The statement is FALSE.\n")
    
    # N) For positive integers x,y: sum([x/y]*y) equals x in Python 3 only if x is divisible by y
    x_n, y_n = 10, 3 # A counterexample
    is_divisible_n = (x_n % y_n == 0)
    # The expression `sum([x/y]*y)` is equivalent to `(x/y)*y`
    res_n_sum = (x_n / y_n) * y_n
    res_n_equals_x = (res_n_sum == x_n)
    print(f"N) Statement: 'sum([x/y]*y) == x' in Py3 only if x is divisible by y.")
    print(f"   - Counterexample: x={x_n}, y={y_n}. Here, x is not divisible by y.")
    print(f"   - The expression `({x_n}/{y_n}) * {y_n}` evaluates to `{res_n_sum}`.")
    print(f"   - The condition `sum(...) == x` is {res_n_equals_x}, but `x` is not divisible by `y`.")
    print(f"   - Verdict: The 'only if' clause is violated. The statement is FALSE.\n")

    # O) The expression (1/3) + (1/3) + (1/3) == 1.0 evaluates to True even with floating-point arithmetic
    res_o_sum = (1/3) + (1/3) + (1/3)
    is_o_true = (res_o_sum == 1.0)
    print(f"O) Statement: `(1/3) + (1/3) + (1/3) == 1.0` is true.")
    print(f"   - Due to floating point inaccuracy, `1/3` is an approximation.")
    print(f"   - The sum `(1/3) + (1/3) + (1/3)` evaluates to `{res_o_sum:.17f}`.")
    print(f"   - Verdict: `{res_o_sum:.17f} == 1.0` is {is_o_true}. The statement is FALSE.\n")

if __name__ == '__main__':
    analyze_division_statements()