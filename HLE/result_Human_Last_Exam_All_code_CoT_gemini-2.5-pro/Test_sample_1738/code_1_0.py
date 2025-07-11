import sys
import io

# Capture original stdout to restore it later if needed, though not strictly necessary here.
original_stdout = sys.stdout
# Create a string stream to capture output for final formatting
output_stream = io.StringIO()
sys.stdout = output_stream

# --- Initial Data ---
a = {1, 2, 3}
b = {3, 4, 5}
c = (a, b)
d = ((1, 2), (3, 4))
true_statements = []

# --- Evaluation of Statements ---

print("--- Evaluating Statement A ---")
# A) The expression c[0] and c[1] always returns the same value as c[1] and c[0]
val1 = c[0] and c[1]
val2 = c[1] and c[0]
print(f"Given a={a}, b={b}, c=(a,b):")
print(f"c[0] and c[1] is ({c[0]}) and ({c[1]}), which evaluates to: {val1}")
print(f"c[1] and c[0] is ({c[1]}) and ({c[0]}), which evaluates to: {val2}")
print(f"Conclusion: The values are not the same. Statement A is False.")

print("\n--- Evaluating Statement B ---")
# B) For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
s1, s2 = {10, 20}, {20, 30}
res_b = (s1 and s2) | (s2 and s1)
print(f"For non-empty sets s1={s1}, s2={s2}:")
print(f"(s1 and s2) evaluates to s2: {s1 and s2}")
print(f"(s2 and s1) evaluates to s1: {s2 and s1}")
print(f"The full expression becomes s2 | s1, which is the union: {res_b}")
print(f"Conclusion: The result {res_b} is not equal to s2 {s2}. Statement B is False.")

print("\n--- Evaluating Statement C ---")
# C) The expression d[0] or d[1] evaluates to (1, 2) regardless of the values in the tuples
d_test_false = ((), (5, 6))
res_c1 = d[0] or d[1]
res_c2 = d_test_false[0] or d_test_false[1]
print(f"For the given d={d}, d[0] or d[1] evaluates to: {res_c1}")
print(f"However, for a different tuple d={d_test_false} where the first element is falsey:")
print(f"The expression evaluates to: {res_c2}")
print(f"Conclusion: The result depends on the truthiness of d[0]. Statement C is False.")

print("\n--- Evaluating Statement D ---")
# D) If s is any non-empty set, s and True returns True while True and s returns s
s = {"apple", "banana"}
res_d1 = s and True
res_d2 = True and s
print(f"For a non-empty set s={s}:")
print(f"s and True evaluates to: {res_d1}")
print(f"True and s evaluates to: {res_d2}")
print("Conclusion: Since s is truthy, 's and True' returns the second operand (True). Since True is truthy, 'True and s' returns the second operand (s). Statement D is True.")
true_statements.append('D')

print("\n--- Evaluating Statement E ---")
# E) The expression {} or [] evaluates to [] while [] or {} evaluates to {}
res_e1 = {} or []
res_e2 = [] or {}
print(f"{{}} or [] evaluates to: {res_e1!r}") # Use !r to distinguish [] from {}
print(f"[] or {{}} evaluates to: {res_e2!r}")
print("Conclusion: In an 'or' chain of falsey values, the last value is returned. Statement E is True.")
true_statements.append('E')

print("\n--- Evaluating Statement F ---")
# F) For any tuple t containing exactly two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
t_test = ({1}, {}) # A truthy and a falsey set
lhs = bool(t_test[0] and t_test[1])
rhs = bool(t_test[0]) and bool(t_test[1])
print(f"Testing with t = {t_test}:")
print(f"LHS: bool(t[0] and t[1]) -> bool({t_test[0] and t_test[1]}) -> {lhs}")
print(f"RHS: bool(t[0]) and bool(t[1]) -> {bool(t_test[0])} and {bool(t_test[1])} -> {rhs}")
print("Conclusion: This holds true for all combinations of truthy/falsey values and is the definition of logical 'and'. Statement F is True.")
true_statements.append('F')

print("\n--- Evaluating Statement G ---")
# G) The expression (a and b) - (b and a) is always empty for any two sets a, b
a_g, b_g = {1, 2}, {2, 3}
res_g = (a_g and b_g) - (b_g and a_g)
print(f"Counterexample with a={a_g}, b={b_g}:")
print(f"(a and b) evaluates to b: {a_g and b_g}")
print(f"(b and a) evaluates to a: {b_g and a_g}")
print(f"The expression becomes b - a, which is: {res_g}")
print(f"Conclusion: The result is not empty. Statement G is False.")

print("\n--- Evaluating Statement H ---")
# H) If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
x_h, y_h = (0, 1), (2, 3) # x[0] is falsey
lhs_h = (x_h and y_h)[0]
rhs_h = x_h[0] and y_h[0]
print(f"Counterexample with x={x_h}, y={y_h}:")
print(f"LHS: (x and y)[0] -> ({x_h and y_h})[0] -> {y_h}[0] -> {lhs_h}")
print(f"RHS: x[0] and y[0] -> {x_h[0]} and {y_h[0]} -> {rhs_h}")
print(f"Conclusion: {lhs_h} is not equal to {rhs_h}. Statement H is False.")

print("\n--- Evaluating Statement I ---")
# I) For any sets p, q, r: (p and q) and r equals p and (q and r)
p, q, r = {1}, {}, {3} # Test case: T, F, T
lhs_i = (p and q) and r
rhs_i = p and (q and r)
print(f"Testing with p={p}, q={q}, r={r}:")
print(f"LHS: (p and q) and r -> ({p and q}) and {r} -> {q} and {r} -> {lhs_i}")
print(f"RHS: p and (q and r) -> {p} and ({q and r}) -> {p} and {q} -> {rhs_i}")
print(f"Conclusion: Both sides evaluate to the first falsey value, or the last value if all are truthy. The `and` operator is associative in this manner. Statement I is True.")
true_statements.append('I')

print("\n--- Evaluating Statement J ---")
# J) The statement all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
try:
    res_j = all(x and y for x, y in zip(a, b))
    print(f"Executing the code `all(x and y for x, y in zip(a, b))`...")
    print(f"The code ran successfully and returned: {res_j}")
    print("Conclusion: No TypeError was raised. Statement J is False.")
except TypeError as e:
    print(f"The code raised a TypeError: {e}")
    print("Conclusion: The statement is True.")

# --- Final Answer ---
true_statements.sort()
final_answer_str = ",".join(true_statements)

# Restore stdout and print the captured output
sys.stdout = original_stdout
print(output_stream.getvalue())
print(f"<<<{final_answer_str}>>>")