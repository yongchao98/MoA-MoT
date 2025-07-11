def solve():
    """
    This function determines for which k in {1, 2, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) vanishes.

    The rational homotopy Lie algebra of a wedge sum is the free product of the
    Lie algebras of the components.
    pi_k(S^4) @ Q is non-zero for k = 4, 7.
    pi_k(CP^2) @ Q is non-zero for k = 2, 5.
    These give non-vanishing groups for X for k in {2, 4, 5, 7}.

    We check the remaining k's for mixed Whitehead products.
    k=1: Vanishes (simply-connected).
    k=3: Vanishes (only product is [beta_2, beta_2], which is 0 as pi_3(CP^2)@Q = 0).
    k=6: Non-zero (e.g., [[beta_2, alpha_4], beta_2]).
    k=8: Non-zero (e.g., [beta_5, alpha_4]).
    k=9: Non-zero (e.g., [[beta_5, alpha_4], beta_2]).

    So, the groups vanish only for k=1 and k=3.
    """
    
    vanishing_k = [1, 3]
    
    # The final answer requires outputting each number in the final list.
    print(','.join(map(str, vanishing_k)))

solve()