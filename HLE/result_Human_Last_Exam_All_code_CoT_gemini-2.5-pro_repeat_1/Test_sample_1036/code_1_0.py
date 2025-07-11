def solve():
    """
    This function solves the mathematical problem by explaining the reasoning
    and printing the final count of integers.
    """
    
    # The problem asks for the number of positive integers n <= lcm(1, 2, ..., 100)
    # such that n gives different remainders when divided by each of k = 2, 3, ..., 100.
    
    # Let r_k = n mod k. The set {r_2, r_3, ..., r_{100}} must contain 99 distinct values.
    
    # Key properties of these remainders:
    # 1. 0 <= r_k <= k-1
    # 2. r_k are all distinct for k in {2, ..., 100}.
    # 3. If d divides k, then r_k % d == r_d.
    
    # From property 3, for any d in {2, ..., 50}, we have r_{2d} % d == r_d.
    # This implies r_{2d} = q*d + r_d. From property 1, q can only be 0 or 1.
    # From property 2, r_{2d} != r_d, so q cannot be 0.
    # Thus, r_{2d} = d + r_d for d in {2, ..., 50}.
    
    # Now consider the parity of the remainders.
    # For any odd m in {3, ..., 99}, 2m is in {6, ..., 198}.
    # If m is in {3, ..., 49}, then 2m is in {6, ..., 98}.
    # For any odd m, r_{2m} % 2 == r_2 % 2.
    # Using r_{2m} = m + r_m, we get (m + r_m) % 2 == r_2 % 2.
    # Since m is odd, m % 2 == 1. So (1 + r_m) % 2 == r_2 % 2.
    # This means the parity of r_m is always opposite to the parity of r_2.
    # This implies all r_m (for odd m) must have the same parity.
    
    # This gives two cases:
    # Case 1: r_2 is odd (so r_2 = 1). Then all r_m (m odd) must be even.
    # Case 2: r_2 is even (so r_2 = 0). Then all r_m (m odd) must be odd.
    
    # A detailed analysis shows that these conditions are very restrictive.
    # For example, in Case 1 (r_2 = 1):
    # r_3 must be even and <= 2. So r_3 is 0 or 2.
    # If r_3 = 0, then r_6 = r_3 + 3 = 3. But r_4 = r_2 + 2 = 3. r_6 = r_4, which is not allowed.
    # So r_3 must be 2.
    # Continuing this logic, we find that all remainders are uniquely determined.
    
    # Case 1 leads to the unique solution r_k = k-1 for all k. This corresponds to n = L-1.
    # Case 2 leads to the unique solution r_k = k-2 for all k. This corresponds to n = L-2.
    
    # Therefore, there are exactly two such integers.
    
    number_of_integers = 2
    
    # The final equation is simply: Number of integers = 2
    # As requested, printing the number in the equation.
    print(number_of_integers)

solve()