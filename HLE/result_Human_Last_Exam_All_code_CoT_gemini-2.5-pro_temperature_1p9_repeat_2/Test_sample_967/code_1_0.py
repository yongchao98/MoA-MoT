def compute_betti_number():
    """
    This script computes the first l^2-Betti number of the fundamental group G
    of the described graph of groups X.

    The computation follows the formula for the first l^2-Betti number of
    the fundamental group of a graph of groups with infinite vertex groups:
    b_1^(2)(G) = sum_{v in V} b_1^(2)(G_v) - sum_{e in E} b_1^(2)(G_e) + b_1(Gamma)
    where Gamma is the underlying graph, G_v are the vertex groups, G_e are
    the edge groups, and b_1(Gamma) is the first Betti number of Gamma.
    """

    # Step 1: Compute the first Betti number b_1(Gamma) of the underlying graph Gamma.
    # Gamma is the line graph of the Petersen graph, L(P).

    # Properties of the Petersen graph (P):
    # It has 10 vertices and 15 edges and is 3-regular.
    petersen_vertices = 10
    petersen_edges = 15

    # The number of vertices in L(P) equals the number of edges in P.
    gamma_vertices = petersen_edges

    # The Petersen graph is 3-regular. The degree of a vertex in L(P) is 
    # (deg(u)-1) + (deg(v)-1) for an edge {u,v}. So L(P) is 4-regular.
    gamma_degree = 4
    # The number of edges in L(P) is (num_vertices * degree) / 2.
    gamma_edges = (gamma_vertices * gamma_degree) // 2

    # The Petersen graph is connected, so its line graph is also connected.
    # The number of connected components b_0(Gamma) is 1.
    gamma_connected_components = 1
    
    # The first Betti number is b_1(Gamma) = |E| - |V| + b_0.
    b1_gamma = gamma_edges - gamma_vertices + gamma_connected_components

    # Step 2: Compute the sum of the first l^2-Betti numbers of the vertex groups.
    # For any g >= 2, M_g is a closed hyperbolic 3-manifold. Its fundamental group
    # pi_1(M_g) has first l^2-Betti number b_1^(2)(pi_1(M_g)) = 0.
    # N_g is a subgroup of pi_1(M_g) of index g.
    # By the formula for finite index subgroups, b_1^(2)(N_g) = g * b_1^(2)(pi_1(M_g)) = 0.
    beta1_Ng = 0

    # The vertex group G_v1 is N_100, so b_1^(2)(G_v1) = 0.
    # The other vertex groups G_vi are free products of N_g groups.
    # For free products of infinite groups, the first l^2-Betti number is additive.
    # b_1^(2)(*_g N_g) = sum_g b_1^(2)(N_g) = sum_g(0) = 0.
    # Therefore, the b_1^(2) of every vertex group is 0.
    sum_beta1_Gv = 0

    # Step 3: Compute the sum of the first l^2-Betti numbers of the edge groups.
    # The edge groups G_e are isomorphic to some N_k.
    # Therefore, b_1^(2)(G_e) = b_1^(2)(N_k) = 0 for some k.
    # Thus, the sum is also 0.
    sum_beta1_Ge = 0

    # Step 4: Combine the results using the formula.
    final_betti_number = sum_beta1_Gv - sum_beta1_Ge + b1_gamma

    print(f"The calculation for the first l^2-Betti number of G is as follows:")
    print(f"b_1^(2)(G) = (Sum of b_1^(2) of vertex groups) - (Sum of b_1^(2) of edge groups) + (First Betti number of the graph)")
    print(f"\n1. First Betti number of the graph Gamma = L(Petersen):")
    print(f"   - Vertices: |V| = {gamma_vertices}")
    print(f"   - Edges: |E| = {gamma_edges}")
    print(f"   - b_1(Gamma) = |E| - |V| + 1 = {gamma_edges} - {gamma_vertices} + 1 = {b1_gamma}")

    print(f"\n2. Sum of b_1^(2) of vertex groups: {sum_beta1_Gv}")
    print(f"   - For each N_g, b_1^(2)(N_g) = 0. Since all G_v are free products of N_g groups (or N_100), their b_1^(2) are also 0.")

    print(f"\n3. Sum of b_1^(2) of edge groups: {sum_beta1_Ge}")
    print(f"   - Each edge group G_e is isomorphic to some N_k, so b_1^(2)(G_e) = 0.")

    print("\nSubstituting these values into the formula:")
    print(f"{sum_beta1_Gv} - {sum_beta1_Ge} + {b1_gamma} = {final_betti_number}")
    
    print(f"\nThe first l^2-Betti number of G is {final_betti_number}.")

compute_betti_number()