def main():
    """
    Computes the first l2-betti number for the given graph of groups.
    """
    
    # 1. Define the structure of the underlying graph: Line graph of the Petersen graph
    # The Petersen graph has 10 vertices (degree 3) and 15 edges.
    # The line graph L(Petersen) has:
    # - Number of vertices = |E(Petersen)| = 15
    # - Number of edges = sum(d(v)*(d(v)-1)/2 for v in V(Petersen)) / 2 --> is not right
    # - Number of edges = 1/2 * sum(deg(e) for e in E(Petersen)) where deg(e) is number of edges incident to e
    #   Each edge in Petersen graph is incident to 4 other edges. Degree of each vertex in L(Petersen) is 4.
    #   So, num_edges = (15 * 4) / 2 = 30
    num_vertices = 15
    num_edges = 30
    
    # 2. Compute the first l2-betti number for each component group N_g
    def get_betti_1_Ng(g):
        """
        Calculates beta_1^{(2)}(N_g).
        N_g is the fundamental group of a g-sheeted covering of M_g.
        M_g is the mapping torus of a pseudo-Anosov map on S_g.
        - M_g is a closed hyperbolic 3-manifold.
        - For a closed hyperbolic 3-manifold M, beta_1^{(2)}(pi_1(M)) = 0.
        - For a finite index subgroup H < K, beta_1^{(2)}(H) = [K:H] * beta_1^{(2)}(K).
        - So, beta_1^{(2)}(N_g) = g * beta_1^{(2)}(pi_1(M_g)) = g * 0 = 0.
        """
        return 0

    # 3. Compute the sum of the first l2-betti numbers of the vertex groups
    betti_v_list = []
    
    # For vertex v_1, the group is N_100
    b_v1 = get_betti_1_Ng(100)
    betti_v_list.append(b_v1)
    
    # For vertices v_i, i>=2, the group is the free product *_{g=2 to i} N_g
    # For infinite groups, beta_1^{(2)} of a free product is the sum of beta_1^{(2)} of factors.
    for i in range(2, num_vertices + 1):
        # The sum of beta_1^{(2)}(N_g) for g from 2 to i
        b_vi = sum(get_betti_1_Ng(g) for g in range(2, i + 1))
        betti_v_list.append(b_vi)
    
    sum_betti_v = sum(betti_v_list)
    
    # 4. Compute the sum of the first l2-betti numbers of the edge groups
    betti_e_list = []
    # An edge group G_e is a "freely indecomposable free factor", meaning G_e is isomorphic to some N_k.
    # The l2-betti number beta_1^{(2)}(N_k) is 0 for any k.
    # Therefore, beta_1^{(2)}(G_e) = 0 for all edges.
    for _ in range(num_edges):
        # We don't need to know which specific N_k it is, the result is always 0.
        b_e = 0
        betti_e_list.append(b_e)
        
    sum_betti_e = sum(betti_e_list)

    # 5. Apply the formula for the l2-betti number of a graph of groups
    # beta_1(G) = sum(beta_1(G_v)) - sum(beta_1(G_e))
    result = sum_betti_v - sum_betti_e

    # 6. Print the detailed equation and the final result
    print("The first l2-betti number of G, beta_1^{(2)}(G), is computed by the formula:")
    print("beta_1^{(2)}(G) = sum(beta_1^{(2)}(G_v)) - sum(beta_1^{(2)}(G_e))")
    print("\nCalculating the sum over vertex groups:")
    v_sum_str = " + ".join(map(str, betti_v_list))
    print(f"sum(beta_1^{(2)}(G_v)) = {v_sum_str} = {sum_betti_v}")

    print("\nCalculating the sum over edge groups:")
    e_sum_str = " + ".join(map(str, betti_e_list))
    print(f"sum(beta_1^{(2)}(G_e)) = {e_sum_str} = {sum_betti_e}")

    print("\nFinal Computation:")
    print(f"beta_1^{(2)}(G) = ({v_sum_str}) - ({e_sum_str})")
    print(f"beta_1^{(2)}(G) = {sum_betti_v} - {sum_betti_e} = {result}")

if __name__ == '__main__':
    main()
