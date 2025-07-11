import sys
import io

def solve_and_print_answer():
    """
    This function analyzes each statement from the problem, determines which are true,
    and prints the sorted letters of the true statements.
    """
    # Given values
    a = {1, 2, 3}
    b = {3, 4, 5}
    c = (a, b)
    d = ((1, 2), (3, 4))

    true_statements = []

    # --- Evaluation of Each Statement ---

    # A: The expression c[0] and c[1] always returns the same value as c[1] and c[0]
    # Test with given values: (a and b) vs (b and a) -> b vs a
    if (c[0] and c[1]) == (c[1] and c[0]):
        # This condition is false for the given non-equal sets a and b.
        true_statements.append('A')

    # B: For any two non-empty sets s1, s2, (s1 and s2) | (s2 and s1) equals s2
    # Test with a counterexample
    s1, s2 = {1}, {2}
    # (s1 and s2) | (s2 and s1) -> s2 | s1 -> {1, 2}, which is not equal to s2
    if ((s1 and s2) | (s2 and s1)) == s2:
        true_statements.append('B')

    # C: The expression d[0] or d[1] evaluates to (1, 2)
    # d[0] is (1, 2) which is truthy. The expression short-circuits to d[0].
    if (d[0] or d[1]) == (1, 2):
        true_statements.append('C')

    # D: If s is any non-empty set, s and True returns True while True and s returns s
    s = {"non-empty"} # A sample truthy set
    if (s and True is True) and (True and s is s):
        true_statements.append('D')

    # E: The expression {} or [] evaluates to [] while [] or {} evaluates to {}
    if (({} or []) == []) and (([] or {}) == {}):
        true_statements.append('E')

    # F: For any tuple t containing two sets, bool(t[0] and t[1]) equals bool(t[0]) and bool(t[1])
    # This is a fundamental property of boolean logic that holds true in Python.
    # Test cases:
    t1 = ({1}, {2})   # truthy, truthy
    t2 = (set(), {2}) # falsey, truthy
    t3 = ({1}, set()) # truthy, falsey
    t4 = (set(), set())# falsey, falsey
    f_is_true = True
    for t_case in [t1, t2, t3, t4]:
        if (bool(t_case[0] and t_case[1])) != (bool(t_case[0]) and bool(t_case[1])):
            f_is_true = False
            break
    if f_is_true:
        true_statements.append('F')

    # G: (a and b) - (b and a) is always empty for any two sets a, b
    # Counterexample from B: a={1}, b={2}. (a and b) - (b and a) -> b - a -> {2}
    g_a, g_b = {1}, {2}
    if len((g_a and g_b) - (g_b and g_a)) == 0:
        true_statements.append('G')

    # H: If x and y are tuples of the same length, (x and y)[0] is always equal to x[0] and y[0]
    # Counterexample: x=(0,1), y=(2,3).
    # (x and y)[0] -> y[0] -> 2.
    # x[0] and y[0] -> 0 and 2 -> 0. Not equal.
    h_x, h_y = (0, 1), (2, 3)
    if (h_x and h_y)[0] == (h_x[0] and h_y[0]):
        true_statements.append('H')

    # I: For any sets p, q, r: (p and q) and r equals p and (q and r)
    # This is true due to the associativity of Python's 'and' operator.
    # We can test with a combination of truthy and falsey sets.
    p, q, r = {1}, set(), {3} # truthy, falsey, truthy
    if ((p and q) and r) == (p and (q and r)):
         true_statements.append('I')

    # J: all(x and y for x, y in zip(a, b)) raises TypeError for sets a, b
    # zip(set, set) is valid. `and` on integers is valid. No TypeError.
    try:
        all(x and y for x, y in zip(a, b))
    except TypeError:
        true_statements.append('J') # This will not be reached.

    # Sort the final list of true statements and print
    true_statements.sort()
    final_answer = ",".join(true_statements)
    print(f"<<<{final_answer}>>>")

# Execute the function to find and print the solution
solve_and_print_answer()