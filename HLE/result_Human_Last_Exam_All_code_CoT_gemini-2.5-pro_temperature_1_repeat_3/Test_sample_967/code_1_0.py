import numpy as np

def compute_l2_betti_number():
    """
    Computes the first l2-Betti number for the given graph of groups.
    """

    # Step 1: Analyze the graph Gamma = L(P)
    # The Petersen graph P has 10 vertices and 15 edges. It is 3-regular.
    # The line graph L(P) has vertices corresponding to edges of P, and edges
    # corresponding to adjacent edges in P.
    num_vertices_P = 10
    num_edges_P = 15

    # Number of vertices in Gamma = L(P) is the number of edges in P.
    V = num_edges_P

    # The Petersen graph is 3-regular. The degree of a vertex in L(P)
    # corresponding to an edge {u,v} in P is deg(u)-1 + deg(v)-1.
    # So, L(P) is (3-1)+(3-1) = 4-regular.
    degree_L_P = 4
    # Number of edges in Gamma is (V * degree) / 2
    E = (V * degree_L_P) // 2

    # First Betti number of Gamma
    betti1_gamma = E - V + 1
    
    print("--- Graph Properties (Gamma = L(P)) ---")
    print(f"Number of vertices |V|: {V}")
    print(f"Number of edges |E|: {E}")
    print(f"First Betti number beta_1(Gamma) = |E| - |V| + 1 = {E} - {V} + 1 = {betti1_gamma}\n")

    # Step 2: Analyze the vertex groups G_v
    # For any g>=2, Ng is an infinite group with beta_1^(2)(Ng) = 0 and beta_0^(2)(Ng) = 0.
    # All vertex groups G_v are infinite, so beta_0^(2)(G_v) = 0 for all v.
    sum_beta0_Gv = 0

    # Calculate beta_1^(2)(G_v) for each vertex v
    # For v_1, G_v1 = N_100, so beta_1^(2)(G_v1) = 0.
    beta1_G_v1 = 0
    
    # For v_i, i >= 2, G_vi is the free product of i-1 infinite groups (N_2, ..., N_i).
    # The formula for the first l2-Betti number of a free product of k infinite groups H_j is:
    # beta_1^(2)(H_1*...*H_k) = sum(beta_1^(2)(H_j)) + k - 1.
    # Here, beta_1^(2)(N_g) = 0 for all g.
    # So, beta_1^(2)(G_vi) = 0 + (i-1 - 1) = i - 2. This works for i >= 2.
    
    sum_beta1_Gv = beta1_G_v1
    # The vertices are enumerated v_1, ..., v_15
    for i in range(2, V + 1):
        beta1_G_vi = i - 2
        sum_beta1_Gv += beta1_G_vi

    vertex_contribution = sum_beta1_Gv - sum_beta0_Gv
    
    print("--- Vertex Group Contributions ---")
    print("For all vertex groups G_v, beta_0^(2)(G_v) = 0 since they are infinite.")
    print(f"Sum of beta_1^(2)(G_v) over all {V} vertices = {sum_beta1_Gv}")
    print(f"Total vertex contribution = sum(beta_1^(2)(G_v) - beta_0^(2)(G_v)) = {sum_beta1_Gv} - {sum_beta0_Gv} = {vertex_contribution}\n")

    # Step 3: Analyze the edge groups G_e
    # Edges connected to v_1: The graph is 4-regular, so there are 4 such edges.
    # G_v1 = N_100 is freely indecomposable. For an edge {v_1, v_j}, the edge group must be a
    # free factor of N_100 and G_vj. This forces the edge group to be trivial.
    # For a trivial group H={1}, |H|=1, so beta_0^(2)(H) = 1/1 = 1.
    num_edges_from_v1 = 4
    beta0_trivial_group = 1
    
    # Other edges: The remaining E - 4 edges connect v_i, v_j with i,j >= 2.
    # The edge group is some N_k, which is infinite. So its beta_0^(2) is 0.
    num_other_edges = E - num_edges_from_v1
    beta0_infinite_group = 0

    edge_contribution = (num_edges_from_v1 * beta0_trivial_group) + (num_other_edges * beta0_infinite_group)
    
    print("--- Edge Group Contributions ---")
    print(f"Number of edges with trivial group: {num_edges_from_v1} (contribution: {beta0_trivial_group} each)")
    print(f"Number of edges with infinite group: {num_other_edges} (contribution: {beta0_infinite_group} each)")
    print(f"Total edge contribution = sum(beta_0^(2)(G_e)) = {num_edges_from_v1} * {beta0_trivial_group} + {num_other_edges} * {beta0_infinite_group} = {edge_contribution}\n")
    
    # Step 4: Final calculation
    # beta_1^(2)(G) = vertex_contribution + edge_contribution + betti1_gamma
    final_betti_number = vertex_contribution + edge_contribution + betti1_gamma
    
    print("--- Final Calculation ---")
    print("beta_1^(2)(G) = (Vertex Contribution) + (Edge Contribution) + beta_1(Gamma)")
    print(f"beta_1^(2)(G) = {vertex_contribution} + {edge_contribution} + {betti1_gamma}")
    print(f"The first l^2-Betti number is: {final_betti_number}")
    
    return final_betti_number

compute_l2_betti_number()