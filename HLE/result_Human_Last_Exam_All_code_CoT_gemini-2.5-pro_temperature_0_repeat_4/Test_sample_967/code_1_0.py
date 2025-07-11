import math

def solve_betti_number():
    """
    This script calculates the first l2-Betti number of the fundamental group G
    of the described graph of groups.
    """

    # Step 1: Analyze the underlying graph (Line Graph of the Petersen Graph)
    # The Petersen graph has 10 vertices and 15 edges. It is 3-regular.
    # The line graph L(P) has vertices corresponding to edges of P, and edges
    # corresponding to adjacent edges in P.
    num_vertices = 15  # Number of edges in Petersen graph
    num_edges = 15     # Sum of C(deg(v), 2) for v in P = 10 * C(3,2) = 10 * 3 = 30. Wait, this is wrong.
                       # The number of edges in the line graph L(G) is (1/2) * sum_{v in V(G)} (deg(v)^2) - |E(G)|
                       # For Petersen graph (3-regular): (1/2) * 10 * (3^2) - 15 = 45 - 15 = 30.
                       # Let me re-check the number of edges of L(P).
                       # Each vertex in P has degree 3, so it corresponds to a triangle in L(P).
                       # There are 10 such vertices, so 10 triangles. Total edges = 10 * 3 = 30.
                       # However, some edges might be shared.
                       # A simpler formula for a k-regular graph is |E(L(G))| = k * |V(G)| / 2.
                       # The line graph of a k-regular graph is (2k-2)-regular.
                       # L(P) is (2*3-2) = 4-regular.
                       # So, num_edges = (num_vertices * degree) / 2 = (15 * 4) / 2 = 30.
    num_edges = 30
    
    # The graph is connected.
    b_0 = 1
    # First Betti number of the graph
    b_1_graph = num_edges - num_vertices + b_0

    # Step 2: Compute the sum of the first l2-Betti numbers of the vertex groups.
    # The vertices are enumerated v_1, ..., v_15.
    # For any g>=2, Ng is a finite-index subgroup of the fundamental group of a
    # closed hyperbolic 3-manifold, so its l2-Betti numbers are all 0.
    # beta_1(N_g) = 0 for all g.
    #
    # The formula for a free product of k infinite groups is:
    # beta_1(G_1 * ... * G_k) = sum(beta_1(G_i)) + k - 1
    
    # For v_1, the group is G_v1 = N_100.
    beta1_Gv1 = 0
    
    # For v_i (i >= 2), the group is G_vi = *_{g=2 to i} N_g.
    # This is a free product of k = i-1 groups.
    # beta_1(G_vi) = sum_{g=2 to i} beta_1(N_g) + (i-1) - 1 = 0 + i - 2.
    
    sum_beta1_Gv = beta1_Gv1
    vertex_contributions = [beta1_Gv1]
    for i in range(2, num_vertices + 1):
        beta1_Gvi = i - 2
        sum_beta1_Gv += beta1_Gvi
        vertex_contributions.append(beta1_Gvi)

    # Step 3: Compute the sum of the first l2-Betti numbers of the edge groups.
    # Edge groups G_e are "freely indecomposable free factors", which are the N_g groups.
    # Since beta_1(N_g) = 0, the contribution from each edge group is 0.
    sum_beta1_Ge = 0

    # Step 4: Apply the main formula for the graph of groups.
    # beta_1(G) = sum_v(beta_1(G_v)) - sum_e(beta_1(G_e)) + b_1(graph)
    final_betti_number = sum_beta1_Gv - sum_beta1_Ge + b_1_graph

    # Step 5: Print the results and the final equation.
    print("### Calculation of the First l2-Betti Number ###\n")
    print("The first l2-Betti number of the fundamental group G is computed using the formula:")
    print("β₁⁽²⁾(G) = Σᵥ β₁⁽²⁾(Gᵥ) - Σₑ β₁⁽²⁾(Gₑ) + b₁(Γ)\n")
    
    print("1. Sum of vertex group contributions (Σᵥ β₁⁽²⁾(Gᵥ)):")
    print(f"   - The vertex group Gᵥ₁ is N₁₀₀, so β₁⁽²⁾(Gᵥ₁) = 0.")
    print(f"   - For vertices vᵢ (i=2,...,15), Gᵥᵢ is a free product of i-1 groups N₉.")
    print(f"   - The formula gives β₁⁽²⁾(Gᵥᵢ) = i - 2.")
    print(f"   - The sum is 0 + (2-2) + (3-2) + ... + (15-2) = Σ_{{j=0}}^{{13}} j.")
    print(f"   - Σᵥ β₁⁽²⁾(Gᵥ) = {sum_beta1_Gv}\n")

    print("2. Sum of edge group contributions (Σₑ β₁⁽²⁾(Gₑ)):")
    print(f"   - Edge groups are factors Nₖ, for which β₁⁽²⁾(Nₖ) = 0.")
    print(f"   - Σₑ β₁⁽²⁾(Gₑ) = {num_edges} * 0 = {sum_beta1_Ge}\n")

    print("3. First Betti number of the underlying graph Γ (b₁(Γ)):")
    print(f"   - Γ is the line graph of the Petersen graph.")
    print(f"   - Number of vertices |V| = {num_vertices}")
    print(f"   - Number of edges |E| = {num_edges}")
    print(f"   - b₁(Γ) = |E| - |V| + 1 = {num_edges} - {num_vertices} + 1 = {b_1_graph}\n")

    print("### Final Equation ###")
    print(f"β₁⁽²⁾(G) = {sum_beta1_Gv} - {sum_beta1_Ge} + {b_1_graph}")
    print(f"β₁⁽²⁾(G) = {final_betti_number}")

solve_betti_number()
<<<107>>>