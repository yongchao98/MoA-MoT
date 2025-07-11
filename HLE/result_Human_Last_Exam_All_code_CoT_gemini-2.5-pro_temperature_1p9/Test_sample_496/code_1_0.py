def solve_adjoint_cohomology():
    """
    Calculates the total rank of A = H_{SO(4)}^*(SO(4) \ X) up to degree 100.
    """
    limit = 100

    # Calculate Sum_{k=0 to 100} rank(H_G^k(G))
    # H_G^*(G; Q) is isomorphic to Q[t_L, t_R] with deg(t_L)=deg(t_R)=2.
    # The rank in degree k=2m is m+1. Rank is 0 for odd degrees.
    sum_g_k = 0
    for k in range(limit + 1):
        if k % 2 == 0:
            m = k // 2
            rank_g_k = m + 1
            sum_g_k += rank_g_k

    # Calculate Sum_{k=0 to 100} rank(H_G^k(X))
    # H_G^*(X; Q) is isomorphic to Q[c] with deg(c)=4.
    # The rank in degree k=4m is 1. Rank is 0 for other degrees.
    sum_x_k = 0
    for k in range(limit + 1):
        if k % 4 == 0:
            rank_x_k = 1
            sum_x_k += rank_x_k
            
    # Calculate Sum_{k=0 to 100} rank(H_G^{k-2}(X))
    # This is equivalent to Sum_{j=-2 to 98} rank(H_G^j(X))
    sum_x_k_minus_2 = 0
    for j in range(-2, limit - 2 + 1):
        if j >= 0 and j % 4 == 0:
            rank_x_j = 1
            sum_x_k_minus_2 += rank_x_j

    total_rank = sum_g_k - sum_x_k + sum_x_k_minus_2
    
    print(f"Based on the Gysin sequence in equivariant cohomology, the total rank of A = H_{SO(4)}^*(SO(4) \\ X) for degrees * <= {limit} is given by:")
    print(f"Total Rank = (Sum of ranks of H_G^k(G)) - (Sum of ranks of H_G^k(X)) + (Sum of ranks of H_G^{{k-2}}(X)) for k from 0 to {limit}.")
    print("\n--- Calculation ---")
    print(f"Sum of ranks of H_G^*(G; Q) up to degree {limit}: {sum_g_k}")
    print(f"Sum of ranks of H_G^*(X; Q) up to degree {limit}: {sum_x_k}")
    print(f"Sum of ranks of H_G^*(X; Q) up to degree {limit-2}: {sum_x_k_minus_2}")
    print("\n--- Final Equation ---")
    print(f"Total Rank = {sum_g_k} - {sum_x_k} + {sum_x_k_minus_2} = {total_rank}")


solve_adjoint_cohomology()
<<<1325>>>