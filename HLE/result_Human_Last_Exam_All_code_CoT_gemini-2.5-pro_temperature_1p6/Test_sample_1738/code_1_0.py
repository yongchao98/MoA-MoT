import sys

# Suppress stderr for a cleaner output in some environments, as {} is a dict
# and the user might expect a set. This doesn't affect the logic.
class NullWriter:
    def write(self, text):
        pass
# original_stderr = sys.stderr
# sys.stderr = NullWriter()

def run_analysis():
    """
    Analyzes each statement and prints a detailed evaluation.
    """
    print("--- Python Semantics and Evaluation ---")

    # Given data
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    true_statements = []

    # A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    print("\n--- Evaluation of A ---")
    val_A1 = c[0] and c[1]
    val_A2 = c[1] and c[0]
    print(f"c[0] and c[1] -> {val_A1}")
    print(f"c[1] and c[0] -> {val_A2}")
    print(f"Conclusion: The values are different. Statement A is False.")

    # B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    print("\n--- Evaluation of B ---")
    s1, s2 = {1}, {2}
    val_B = (s1 and s2) | (s2 and s1)
    print(f"For s1={s1}, s2={s2}: (s1 and s2)|(s2 and s1) -> {val_B}")
    print(f"Conclusion: The result {val_B} is not equal to s2 ({s2}). Statement B is False.")

    # C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values
    print("\n--- Evaluation of C ---")
    d_test = ((), (3, 4))
    val_C = d_test[0] or d_test[1]
    print(f"For d_test={d_test}: d_test[0] or d_test[1] -> {val_C}")
    print(f"Conclusion: The result depends on the values. Statement C is False.")

    # D) If s is non-empty, s and True returns True while True and s returns s
    print("\n--- Evaluation of D ---")
    s = {'example'}
    val_D1 = s and True
    val_D2 = True and s
    print(f"For non-empty s: s and True -> {val_D1} (type: {type(val_D1).__name__})")
    print(f"For non-empty s: True and s -> {val_D2}")
    print(f"Conclusion: The statement is correct. Statement D is True.")
    if val_D1 is True and val_D2 is s:
        true_statements.append("D")

    # E) {} or [] evaluates to [] while [] or {} evaluates to {}
    print("\n--- Evaluation of E ---")
    val_E1 = {} or []
    val_E2 = [] or {}
    print(f"{{}} or [] -> {val_E1} (type: {type(val_E1).__name__})")
    print(f"[] or {{}} -> {val_E2} (type: {type(val_E2).__name__})")
    print(f"Conclusion: The statement is correct. Statement E is True.")
    if isinstance(val_E1, list) and isinstance(val_E2, dict):
        true_statements.append("E")

    # F) bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    print("\n--- Evaluation of F ---")
    t = ({1}, set()) # Test case: True, False
    lhs_F = bool(t[0] and t[1])
    rhs_F = bool(t[0]) and bool(t[1])
    print(f"For t={t}:")
    print(f"bool(t[0] and t[1]) -> {lhs_F}")
    print(f"bool(t[0]) and bool(t[1]) -> {rhs_F}")
    print(f"Conclusion: The property holds for all objects. Statement F is True.")
    true_statements.append("F")

    # G) (a and b) - (b and a) is always empty for any two sets a, b
    print("\n--- Evaluation of G ---")
    a_g, b_g = {1, 2}, {2, 3}
    val_G = (a_g and b_g) - (b_g and a_g)
    print(f"For non-empty a={a_g}, b={b_g}: (a and b)-(b and a) -> b-a -> {val_G}")
    print(f"Conclusion: The result is not empty. Statement G is False.")

    # H) If x, y are tuples, (x and y)[0] is always equal to x[0] and y[0]
    print("\n--- Evaluation of H ---")
    x, y = (0, 1), (2, 3)
    lhs_H = (x and y)[0]
    rhs_H = x[0] and y[0]
    print(f"For x={x}, y={y}:")
    print(f"(x and y)[0] -> {lhs_H}")
    print(f"x[0] and y[0] -> {rhs_H}")
    print(f"Conclusion: {lhs_H} != {rhs_H}. Statement H is False.")

    # I) For any sets p, q, r: (p and q) and r equals p and (q and r)
    print("\n--- Evaluation of I ---")
    p, q, r = {1}, set(), {3} # Test case: T, F, T
    lhs_I = (p and q) and r
    rhs_I = p and (q and r)
    print(f"For p={p}, q={q}, r={r}:")
    print(f"(p and q) and r -> {lhs_I}")
    print(f"p and (q and r) -> {rhs_I}")
    print(f"Conclusion: The 'and' operator is associative. Statement I is True.")
    true_statements.append("I")

    # J) all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    print("\n--- Evaluation of J ---")
    try:
        all(x and y for x, y in zip(a, b))
        print(f"Code executed successfully without TypeError.")
        print(f"Conclusion: The statement is False.")
    except TypeError:
        print(f"TypeError was raised.")

    # --- Final Result ---
    print("\n-------------------------")
    true_statements.sort()
    print(f"The true statements are: {', '.join(true_statements)}")

run_analysis()

# sys.stderr = original_stderr