def solve_betti_number():
    """
    Computes the first l2-betti number of the fundamental group G.
    """

    # Step 1: Compute the graph-theoretic term of the formula.
    # The graph X is the line graph of the Petersen graph.
    # The Petersen graph has 10 vertices and 15 edges.
    # The number of vertices in X is the number of edges in the Petersen graph.
    num_vertices_X = 15
    # The number of edges in X (a line graph of a 3-regular graph on 10 vertices) is 30.
    num_edges_X = 30
    
    # The formula component from the graph structure is |E(X)| - |V(X)|.
    # This equals beta_1(X) - b_0(X).
    graph_term = num_edges_X - num_vertices_X

    # Step 2: Compute the sum of the first l2-Betti numbers of the edge groups.
    # Each edge group G_e is some N_g.
    # For a hyperbolic 3-manifold M_g that fibers over the circle, beta_1^(2)(pi_1(M_g)) = 0.
    # N_g is a finite index subgroup, so beta_1^(2)(N_g) is also 0.
    sum_beta_1_Ge = 0

    # Step 3: Compute the sum of the first l2-Betti numbers of the vertex groups.
    # For v_1, the group is N_100, so its beta_1^(2) is 0.
    beta_1_Gv1 = 0

    # For v_i (i from 2 to 15), the group G_{v_i} is a free product of i-1 infinite groups N_g.
    # The formula for the first l2-betti number of a free product of k infinite groups
    # is the sum of their individual l2-betti numbers plus k-1.
    # Since beta_1^(2)(N_g) = 0, for G_{v_i} this becomes 0 + (i-1) - 1 = i - 2.
    # We sum this from i = 2 to 15.
    num_total_vertices = 15
    sum_beta_1_Gvi_for_i_ge_2 = sum(i - 2 for i in range(2, num_total_vertices + 1))

    # The total sum for vertex groups:
    sum_beta_1_Gv = beta_1_Gv1 + sum_beta_1_Gvi_for_i_ge_2

    # Step 4: Combine the terms using the formula:
    # beta_1^(2)(G) = sum_v(beta_1^(2)(G_v)) - sum_e(beta_1^(2)(G_e)) + (|E(X)| - |V(X)|)
    final_betti_number = sum_beta_1_Gv - sum_beta_1_Ge + graph_term

    # Print the final equation with the computed values.
    print(f"The first l2-betti number is the result of the following calculation:")
    print(f"{sum_beta_1_Gv} - {sum_beta_1_Ge} + {graph_term} = {final_betti_number}")

solve_betti_number()