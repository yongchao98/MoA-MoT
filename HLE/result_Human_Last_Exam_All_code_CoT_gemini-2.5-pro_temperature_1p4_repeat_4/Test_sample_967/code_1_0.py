def compute_l2_betti_number():
    """
    This script computes the first l^2-Betti number for the fundamental group
    of the described graph of groups.
    """
    # Step 1: Define properties of the Petersen graph
    petersen_vertices = 10
    petersen_edges = 15
    petersen_degree = 3  # The Petersen graph is 3-regular

    # Step 2: Calculate properties of the line graph X
    # The vertices of the line graph correspond to the edges of the original graph.
    X_vertices = petersen_edges
    # The number of edges in the line graph L(G) of a k-regular graph G is |E(G)| * (k - 1).
    X_edges = petersen_edges * (petersen_degree - 1)
    # The line graph of a connected graph (that is not a path) is connected.
    X_b0 = 1
    # First Betti number of the graph X
    X_beta1 = X_edges - X_vertices + X_b0

    # Step 3: Determine the contribution from vertex and edge groups.
    # As explained in the plan, the first l^2-Betti numbers of all vertex and edge
    # groups are 0.
    # beta_1^(2)(N_g) = 0 for all g.
    # G_v1 = N_100, so beta_1^(2)(G_v1) = 0.
    # G_vi = ast_{g=2 to i} N_g, so beta_1^(2)(G_vi) = sum(beta_1^(2)(N_g)) = 0.
    # G_e is some N_k, so beta_1^(2)(G_e) = 0.
    sum_beta1_Gv = 0
    sum_beta1_Ge = 0

    # Step 4: Apply the formula for the first l^2-Betti number of a graph of groups
    # Formula: result = (sum_beta1_Gv - sum_beta1_Ge) + X_beta1 - X_b0
    result = (sum_beta1_Gv - sum_beta1_Ge) + X_beta1 - X_b0
    
    # Output the result and the equation
    print("The first l^2-Betti number, beta_1^{(2)}(G), is calculated as follows:")
    print("beta_1^{(2)}(G) = (sum[beta_1^{(2)}(G_v)] - sum[beta_1^{(2)}(G_e)]) + beta_1(X) - b_0(X)")
    print(f"The components of the formula are:")
    print(f"  - First Betti number of the graph, beta_1(X) = {X_edges} - {X_vertices} + {X_b0} = {X_beta1}")
    print(f"  - Number of connected components of the graph, b_0(X) = {X_b0}")
    print(f"  - Sum of first l^2-Betti numbers of vertex groups, sum[beta_1^{(2)}(G_v)] = {sum_beta1_Gv}")
    print(f"  - Sum of first l^2-Betti numbers of edge groups, sum[beta_1^{(2)}(G_e)] = {sum_beta1_Ge}")
    print("\nSubstituting these values into the formula:")
    print(f"beta_1^{(2)}(G) = ({sum_beta1_Gv} - {sum_beta1_Ge}) + {X_beta1} - {X_b0}")
    print(f"beta_1^{(2)}(G) = {sum_beta1_Gv - sum_beta1_Ge} + {X_beta1 - X_b0}")
    print(f"beta_1^{(2)}(G) = {result}")

compute_l2_betti_number()