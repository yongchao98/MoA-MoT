def solve():
    """
    This function explains the logic and calculates the maximum number of real coins
    that can be guaranteed to be identified.
    """
    
    # Let n be the size of groups A, B, C, D.
    # Let E be the size of the remaining group.
    # Total coins = 4 * n + E = 1000
    
    # Strategy:
    # W1: (A+B) vs (C+D)
    # W2: (A+C) vs (B+D)
    
    # We analyze the number of guaranteed real coins for each of the 9 outcomes.
    # The number of guaranteed coins for the strategy is the minimum of these 9 values.
    # Let's denote the size of group X as |X|. |A|=|B|=|C|=|D|=n. |E|=1000-4n.
    
    # Outcome Analysis (simplified):
    # 1. (=,=): min(|A|+|D|, |B|+|C|, |E|) = min(2n, |E|)
    # 2. (<,<): C and D are guaranteed real. |C|+|D| = 2n
    # 3. (<,=): D and E are guaranteed real. |D|+|E| = n + |E|
    
    # The guaranteed number of coins is min( min(2n, 1000-4n), 2n, n + (1000-4n) )
    # which simplifies to min(2n, 1000-4n, 1000-3n)
    
    # We want to find n that maximizes this minimum value.
    # We are looking for max(min(2n, 1000-4n)).
    # The other term 1000-3n will be larger for the optimal n.
    # To maximize min(X,Y), we set X=Y.
    # 2n = 1000 - 4n
    # 6n = 1000
    # n = 1000 / 6 = 166.66...
    
    # Since n must be an integer, we test n=166 and n=167.
    n1 = 166
    guaranteed1 = min(2 * n1, 1000 - 4 * n1)
    
    n2 = 167
    guaranteed2 = min(2 * n2, 1000 - 4 * n2)
    
    # The maximum of these two is our answer.
    max_guaranteed = max(guaranteed1, guaranteed2)
    
    # Let's use n=166 to present the final numbers.
    n = 166
    num_a = n
    num_b = n
    num_c = n
    num_d = n
    num_e = 1000 - 4 * n
    
    print("The strategy involves dividing 1000 coins into 5 groups (A, B, C, D, E).")
    print(f"Group sizes: |A|={num_a}, |B|={num_b}, |C|={num_c}, |D|={num_d}, |E|={num_e}")
    print("\nThe weighings are:")
    print("Weighing 1: (A+B) vs (C+D)")
    print("Weighing 2: (A+C) vs (B+D)")
    print("\nThis strategy guarantees a minimum number of identified real coins, which is maximized when n=166 or 167.")
    print(f"For n={n}, the number of coins we can guarantee to identify is min(2 * {n}, 1000 - 4 * {n}) = min({2*n}, {1000-4*n}) = {max_guaranteed}.")
    print(f"\nFinal Answer: The maximum number of real coins you can guarantee to identify is {max_guaranteed}.")

solve()