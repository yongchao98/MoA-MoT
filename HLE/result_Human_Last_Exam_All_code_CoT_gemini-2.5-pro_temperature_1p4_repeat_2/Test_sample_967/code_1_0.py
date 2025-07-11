import math

def compute_betti_number():
    """
    Computes the first l2-Betti number for the described group G.
    """
    # Step 1: Analyze the underlying graph Gamma = L(P)
    # The Petersen graph P has 10 vertices and is 3-regular.
    petersen_vertices = 10
    petersen_degree = 3
    
    # The number of vertices in the line graph L(P) is the number of edges in P.
    gamma_vertices = (petersen_vertices * petersen_degree) // 2
    
    # The number of edges in the line graph L(P) for a k-regular graph on n vertices
    # is n * k * (k-1) / 2.
    gamma_edges = (petersen_vertices * petersen_degree * (petersen_degree - 1)) // 2
    
    # The Petersen graph is connected, so its line graph is also connected.
    gamma_b0 = 1
    
    # The first Betti number of Gamma.
    gamma_b1 = gamma_edges - gamma_vertices + gamma_b0
    
    print("Step 1: Compute the first Betti number of the underlying graph Gamma.")
    print(f"The graph Gamma is the line graph of the Petersen graph.")
    print(f"Number of vertices in Gamma: |V| = {gamma_vertices}")
    print(f"Number of edges in Gamma: |E| = {gamma_edges}")
    print(f"First Betti number of Gamma: b_1(Gamma) = |E| - |V| + 1 = {gamma_edges} - {gamma_vertices} + 1 = {gamma_b1}")
    
    # Step 2: Compute the sum of the l2-Betti numbers of the vertex groups.
    sum_beta1_Gv = 0
    # For vertex v_1, the group is N_100. beta_1(N_100) = 0.
    beta1_Gv1 = 0
    sum_beta1_Gv += beta1_Gv1
    
    # For vertices v_i, i = 2, ..., 15.
    # The group G_vi is the free product of i-1 infinite groups (N_2, ..., N_i).
    # Each N_g has beta_1(N_g) = 0.
    # The formula is sum(beta_1) + (num_factors - 1).
    # So, beta_1(G_vi) = 0 + ((i-1) - 1) = i - 2 for i > 1.
    for i in range(2, gamma_vertices + 1):
        num_factors = i - 1
        if num_factors > 0:
            sum_beta1_Gv += (num_factors - 1)
            
    print("\nStep 2: Compute the sum of the first l2-Betti numbers for the vertex groups.")
    print("For v_1, beta_1(G_v1) = 0.")
    print("For v_i (i>=2), G_vi is a free product of i-1 groups N_g, so beta_1(G_vi) = i - 2.")
    print(f"The sum is Sum_v beta_1(G_v) = 0 + Sum_{i=2 to 15} (i-2) = {sum_beta1_Gv}")

    # Step 3: Compute the sum of the l2-Betti numbers of the edge groups.
    # Each edge group G_e is some N_g, and beta_1(N_g) = 0.
    sum_beta1_Ge = 0
    
    print("\nStep 3: Compute the sum of the first l2-Betti numbers for the edge groups.")
    print("Each edge group is some N_g, with beta_1(N_g) = 0.")
    print(f"The sum is Sum_e beta_1(G_e) = {gamma_edges} * 0 = {sum_beta1_Ge}")

    # Step 4: Apply the main formula.
    betti1_G = gamma_b1 + sum_beta1_Gv - sum_beta1_Ge
    
    print("\nStep 4: Combine results using the formula for a graph of groups.")
    print("beta_1(G) = b_1(Gamma) + Sum_v beta_1(G_v) - Sum_e beta_1(G_e)")
    print(f"beta_1(G) = {gamma_b1} + {sum_beta1_Gv} - {sum_beta1_Ge}")
    print(f"beta_1(G) = {betti1_G}")

compute_betti_number()