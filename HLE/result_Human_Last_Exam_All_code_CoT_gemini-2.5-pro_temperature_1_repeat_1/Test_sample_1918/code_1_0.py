def max_complementary_choi_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of a complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space of the channel.
        n (int): Dimension of the output Hilbert space of the channel.
        r (int): Rank of the Choi matrix of the channel.

    Returns:
        int: The maximal rank of the Choi matrix of the complementary channel.
    """
    # The maximal rank is constrained by three quantities:
    # 1. The dimension of the channel's output space, n.
    # 2. The sum of the input dimension and the Choi rank, d + r.
    # 3. The product of the input dimension and the Choi rank, d * r.
    
    # Calculate the terms for the minimum function
    d_plus_r = d + r
    d_times_r = d * r
    
    # The maximal rank is the minimum of these three values.
    max_rank = min(n, d_plus_r, d_times_r)
    
    print(f"Given parameters: d = {d}, n = {n}, r = {r}")
    print("The maximal rank of the complementary channel's Choi matrix is given by the formula: min(n, d + r, d * r)")
    print("\nCalculating each term in the formula:")
    print(f"n = {n}")
    print(f"d + r = {d} + {r} = {d_plus_r}")
    print(f"d * r = {d} * {r} = {d_times_r}")
    print("\nThe final result is:")
    print(f"max_rank = min({n}, {d_plus_r}, {d_times_r}) = {max_rank}")
    
    return max_rank

if __name__ == '__main__':
    # Example 1:
    print("--- Example 1 ---")
    d1, n1, r1 = 4, 10, 3
    max_complementary_choi_rank(d1, n1, r1)
    print("\n" + "="*20 + "\n")

    # Example 2:
    print("--- Example 2 ---")
    d2, n2, r2 = 2, 5, 2
    max_complementary_choi_rank(d2, n2, r2)
    print("\n" + "="*20 + "\n")
    
    # Example 3 (where d*r is the minimum):
    print("--- Example 3 ---")
    d3, n3, r3 = 5, 20, 1
    max_complementary_choi_rank(d3, n3, r3)