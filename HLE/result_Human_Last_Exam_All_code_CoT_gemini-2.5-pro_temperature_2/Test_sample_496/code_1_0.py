def solve_total_rank():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4)\X
    for degrees up to 100.

    The method relies on the formula for the rank a_k in each degree k:
    a_k = c_k + d_{k-2}
    where:
    c_k is the rank of H_SO(4)^k(SO(4)), which is 1 for k=0 and 0 otherwise.
    d_j is the rank of H_SO(4)^j(X), which is 1 for j=4m (m>=0) and 0 otherwise.

    The total rank is the sum of a_k for k from 0 to 100.
    """

    # Contribution from the c_k terms (H_SO(4)^k(SO(4)))
    # c_k is non-zero only for k=0, where its rank is 1.
    # For degrees k from 0 to 100, this term contributes a total rank of 1.
    rank_from_G = 1

    # Contribution from the d_{k-2} terms (H_SO(4)^j(X))
    # We need to count how many k in [0, 100] satisfy the condition that
    # k-2 is a non-negative multiple of 4.
    # Let k-2 = 4m, for m >= 0.
    # k = 4m + 2
    # We need 0 <= 4m + 2 <= 100.
    # -2 <= 4m <= 98
    # -0.5 <= m <= 24.5
    # So integer m can be 0, 1, ..., 24.
    # The number of such values of m is 24 - 0 + 1 = 25.
    rank_from_X_complement = 0
    m_max = 24
    rank_from_X_complement = m_max - 0 + 1

    # Total rank is the sum of the two contributions.
    total_rank = rank_from_G + rank_from_X_complement
    
    print("The total rank is the sum of two parts:")
    print(f"1. The rank from the cohomology of SO(4) itself: {rank_from_G}")
    print(f"2. The rank from the cohomology of the submanifold X (shifted): {rank_from_X_complement}")
    print("The final calculation is:")
    print(f"{rank_from_G} + {rank_from_X_complement} = {total_rank}")

solve_total_rank()