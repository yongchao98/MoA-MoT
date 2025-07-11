def solve_betti_number():
    """
    Computes the first l2-betti number for the specified graph of groups G.

    The overall formula is:
    beta_1(G) = sum_{v in V} beta_1(G_v) - sum_{e in E} beta_1(G_e) + b_1(Gamma)

    where Gamma is the underlying graph L(P).
    """

    # 1. Properties of the underlying graph Gamma = L(P)
    # The Petersen graph P has 15 edges. The line graph L(P) has vertices
    # corresponding to these edges.
    num_vertices = 15
    # The Petersen graph is 3-regular, so its line graph is 4-regular.
    # The number of edges is (num_vertices * degree) / 2.
    num_edges = (num_vertices * 4) // 2
    # The first Betti number of a connected graph is |E| - |V| + 1.
    betti1_gamma = num_edges - num_vertices + 1

    # 2. Sum of l2-Betti numbers of vertex groups
    # For v_1, G_1 = N_100, beta_1(G_1) = 0.
    # For v_i (i>=2), G_i is a free product of i-1 groups N_g,
    # each with beta_1(N_g) = 0.
    # beta_1(G_i) = sum(beta_1(N_g)) + (num_factors - 1) = 0 + (i-1) - 1 = i - 2.
    sum_beta1_Gv = 0
    # The loop calculates sum_{i=2 to 15} (i-2), which is 0+1+...+13
    for i in range(2, num_vertices + 1):
        sum_beta1_Gv += (i - 2)
    # The contribution from G_1 (i=1) is 0, so the sum is correct.
    
    # 3. Sum of l2-Betti numbers of edge groups
    # All edge groups are of the form N_k for some k.
    # beta_1(N_k) = 0 for all k.
    sum_beta1_Ge = 0

    # 4. Final Calculation
    result = sum_beta1_Gv - sum_beta1_Ge + betti1_gamma

    # Output the result as an equation
    print(f"{sum_beta1_Gv} - {sum_beta1_Ge} + {betti1_gamma} = {result}")

solve_betti_number()