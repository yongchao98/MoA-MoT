def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of the complement
    of a 3-dimensional closed invariant submanifold X in SO(4), for degrees up to 100.
    
    The calculation is based on the following steps:
    1. The Poincare series for the equivariant cohomology of SO(4) under the adjoint action is
       P_S(t) = (1 + t^3) / (1 - t^4)^2.
       This means that the rank of the cohomology in degree 4j is j+1,
       and in degree 4j+3 is j+1.
    2. The Gysin long exact sequence for the pair (SO(4), SO(4)\X) implies that,
       under a simplifying assumption that the connecting homomorphisms are zero,
       rank(A^k) = rank(S^k) + rank(H_G^{5-k}(X)), where A is the ring we are interested in,
       S is H_G^*(SO(4)), and d=codim(X)=3.
    3. The rank of H_G^j(X) is 1 for all j >= 0, based on a model of X as an S^2-bundle
       over S^1.
    4. We sum these ranks up to degree 100.
    """
    
    limit_degree = 100
    
    # Part 1: Sum of ranks for S^* = H_{SO(4)}^*(SO(4))
    
    # Ranks for degrees of the form 4j
    # 4j <= 100  => j <= 25
    max_j_1 = limit_degree // 4
    # The ranks are j+1. Sum is sum_{j=0}^{max_j_1} (j+1)
    sum_ranks_S1 = (max_j_1 + 1) * (max_j_1 + 2) // 2
    
    # Ranks for degrees of the form 4j+3
    # 4j+3 <= 100 => 4j <= 97 => j <= 24.25
    max_j_2 = (limit_degree - 3) // 4
    # The ranks are j+1. Sum is sum_{j=0}^{max_j_2} (j+1)
    sum_ranks_S2 = (max_j_2 + 1) * (max_j_2 + 2) // 2
    
    total_rank_S = sum_ranks_S1 + sum_ranks_S2
    
    # Part 2: Sum of ranks for the contribution from X
    # We need to sum rank(H_G^{5-k}(X)) for k from 0 to 100.
    # rank(H_G^j(X)) = 1 for j >= 0, and 0 otherwise.
    # So we sum 1 for all k such that 5-k >= 0, i.e., k <= 5.
    sum_ranks_X_contrib = 0
    for k in range(limit_degree + 1):
        if 5 - k >= 0:
            sum_ranks_X_contrib += 1
            
    # Total rank is the sum of the two parts
    total_rank = total_rank_S + sum_ranks_X_contrib
    
    # Print the equation
    print(f"Total Rank = (Sum of ranks for S^4j) + (Sum of ranks for S^4j+3) + (Sum of ranks from X)")
    print(f"Total Rank = ({sum_ranks_S1}) + ({sum_ranks_S2}) + ({sum_ranks_X_contrib})")
    print(f"Total Rank = {total_rank_S} + {sum_ranks_X_contrib}")
    print(f"Total Rank = {total_rank}")

solve()