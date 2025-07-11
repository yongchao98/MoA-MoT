def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring
    of SO(4) \ X in degrees up to 100.
    """
    total_rank = 0
    max_degree = 100

    print("Calculating the total rank by summing the rank at each even degree up to 100.")
    print("The rank at odd degrees is 0.")
    print("-" * 30)

    # Loop through all degrees from 0 to 100
    for k in range(max_degree + 1):
        # Ranks are non-zero only for even degrees
        if k % 2 == 0:
            m = k // 2

            # Rank of H_G^k(G)
            # For k=2m, this is m+1
            rank_H_G_k_G = m + 1

            # Rank of H_G^{k-2}(X)
            # This is non-zero only if k-2 is even and non-negative
            if k >= 2:
                rank_H_G_k_minus_2_X = 1
            else:
                rank_H_G_k_minus_2_X = 0

            # The rank of A^k is the sum of the two ranks
            rank_A_k = rank_H_G_k_G + rank_H_G_k_minus_2_X
            
            print(f"Rank at degree {k}: {rank_A_k}")
            total_rank += rank_A_k

    print("-" * 30)
    print(f"The final sum of the ranks is the total rank.")
    print(f"Total rank for degrees <= {max_degree}: {total_rank}")

solve()
<<<1376>>>