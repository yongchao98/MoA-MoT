import math

def compute_l2_betti_number():
    """
    Computes the first l2-Betti number of the fundamental group G.
    """
    # Step 1: Compute properties of the underlying graph Gamma.
    # Gamma is the line graph of the Petersen graph.
    # The Petersen graph has 10 vertices and 15 edges.
    num_vertices_gamma = 15
    # The Petersen graph is 3-regular, so its line graph is 4-regular.
    # Number of edges = (num_vertices * degree) / 2
    num_edges_gamma = (num_vertices_gamma * 4) // 2

    # The graph is connected, so its first Betti number is |E| - |V| + 1.
    b1_gamma = num_edges_gamma - num_vertices_gamma + 1

    # Step 2: Compute the sum over vertex groups.
    # The formula term is sum_{v in V} (cost(G_v) - 1).
    # The vertices are enumerated v_1, ..., v_15.
    
    # For v_1, G_v1 = N_100. cost(N_100) = 1.
    cost_G_v1 = 1
    term_v1 = cost_G_v1 - 1
    
    sum_cost_v_minus_1 = term_v1
    
    # For v_i (i=2..15), G_vi = *_{g=2 to i} N_g.
    # cost(G_vi) = sum_{g=2 to i} cost(N_g) = sum_{g=2 to i} 1 = i - 1.
    # The term is cost(G_vi) - 1 = (i - 1) - 1 = i - 2.
    sum_other_v_terms = 0
    for i in range(2, num_vertices_gamma + 1):
        sum_other_v_terms += (i - 2)
        
    sum_cost_v_minus_1 += sum_other_v_terms

    # Step 3: Compute the sum over edge groups.
    # For any edge e, G_e is some N_k, so cost(G_e) = 1.
    # The term cost(G_e) - 1 is always 0.
    sum_cost_e_minus_1 = 0

    # Step 4: Apply Gaboriau's formula for cost.
    # cost(G) - 1 = sum_v(cost(G_v)-1) - sum_e(cost(G_e)-1) + b1(Gamma)
    cost_G_minus_1 = sum_cost_v_minus_1 - sum_cost_e_minus_1 + b1_gamma

    # The fundamental group G is infinite, so beta_0^(2)(G) = 0.
    # cost(G) = beta_1^(2)(G) - beta_0^(2)(G) + 1 = beta_1^(2)(G) + 1.
    # Therefore, beta_1^(2)(G) = cost(G) - 1.
    beta_1_l2 = cost_G_minus_1

    print("The first l2-Betti number is computed using Gaboriau's formula for the cost of a group:")
    print("cost(G) - 1 = sum_v(cost(G_v) - 1) - sum_e(cost(G_e) - 1) + b_1(Gamma)")
    print("This can be translated to l2-Betti numbers, yielding:")
    print(f"beta_1^(2)(G) = sum_v(cost(G_v) - 1) - sum_e(cost(G_e) - 1) + b_1(Gamma)")
    print("-" * 20)
    print(f"1. Graph properties:")
    print(f"   - Number of vertices |V| in Gamma = {num_vertices_gamma}")
    print(f"   - Number of edges |E| in Gamma = {num_edges_gamma}")
    print(f"   - First Betti number b_1(Gamma) = |E| - |V| + 1 = {b1_gamma}")
    print("-" * 20)
    print(f"2. Sum over vertex groups:")
    print(f"   - The term for v_1 is {term_v1}.")
    print(f"   - The sum of terms for v_2 to v_15 is {sum_other_v_terms}.")
    print(f"   - Total sum_v(cost(G_v) - 1) = {sum_cost_v_minus_1}")
    print("-" * 20)
    print(f"3. Sum over edge groups:")
    print(f"   - For each edge e, cost(G_e) = 1, so cost(G_e) - 1 = 0.")
    print(f"   - Total sum_e(cost(G_e) - 1) = {sum_cost_e_minus_1}")
    print("-" * 20)
    print(f"4. Final calculation:")
    print(f"beta_1^(2)(G) = {sum_cost_v_minus_1} - {sum_cost_e_minus_1} + {b1_gamma}")
    print(f"beta_1^(2)(G) = {beta_1_l2}")

compute_l2_betti_number()