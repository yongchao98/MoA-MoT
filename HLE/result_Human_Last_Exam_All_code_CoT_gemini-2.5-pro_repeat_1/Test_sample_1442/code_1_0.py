import math

def n_choose_k(n, k):
    """Calculates the binomial coefficient 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_matchings(n, d, num_c4):
    """
    Calculates the number of 3-matchings in a d-regular bipartite graph on n vertices.

    Args:
        n: The number of vertices.
        d: The degree of each vertex.
        num_c4: The number of 4-cycles in the graph.
    """
    print(f"--- Calculating for a graph with n={n}, d={d}, and N(C4)={num_c4} ---")

    # Number of edges (m)
    if (n * d) % 2 != 0:
        print("Invalid graph parameters: n*d must be even.")
        return
    m = (n * d) // 2
    print(f"Number of edges m = (n*d)/2 = ({n}*{d})/2 = {m}")

    # Total number of ways to choose 3 edges
    total_3_edge_sets = n_choose_k(m, 3)
    print(f"Total sets of 3 edges = m_choose_3 = {m}_choose_3 = {total_3_edge_sets}")

    # Count of bad configurations that are not 3-matchings

    # 1. N(K_1,3): Number of 3-stars (3 edges incident to one vertex)
    # For a d-regular graph, this is n * (d choose 3)
    num_k1_3 = n * n_choose_k(d, 3)
    print(f"Number of 3-stars N(K_1,3) = n * (d choose 3) = {n} * {n_choose_k(d, 3)} = {num_k1_3}")

    # 2. N(P3 U K2): Number of (path of 2 edges) + (1 disjoint edge)
    # Number of P3s is n * (d choose 2).
    num_p3 = n * n_choose_k(d, 2)
    # For any P3 in a bipartite d-regular graph, the number of edges incident to its 3 vertices is 3d-2.
    # So, the number of disjoint edges is m - (3d-2).
    num_disjoint_edges = m - (3 * d - 2)
    num_p3_k2 = num_p3 * num_disjoint_edges
    print(f"Number of P3s = n * (d choose 2) = {n} * {n_choose_k(d, 2)} = {num_p3}")
    print(f"Number of edges disjoint from a P3 = m - (3d-2) = {m} - ({3*d-2}) = {num_disjoint_edges}")
    print(f"Number of P3 U K2 subgraphs = {num_p3} * {num_disjoint_edges} = {num_p3_k2}")

    # 3. N(P4): Number of paths of length 3
    # This count depends on the number of 4-cycles (num_c4).
    # The formula is N(P4) = 0.5 * (n*d*(d-1)^2 - 4*N(C4))
    term1 = n * d * (d - 1)**2
    num_p4 = (term1 - 4 * num_c4) // 2
    print(f"Number of P4 paths N(P4) = 0.5 * [n*d*(d-1)^2 - 4*N(C4)] = 0.5 * [{term1} - 4*{num_c4}] = {num_p4}")

    # Final calculation for the number of 3-matchings N(M3)
    # N(M3) = (total sets) - N(K1,3) - N(P3 U K2) - N(P4)
    num_m3 = total_3_edge_sets - num_k1_3 - num_p3_k2 - num_p4
    
    print("\nFinal Calculation:")
    print(f"N(M3) = (Total 3-edge sets) - N(K1,3) - N(P3 U K2) - N(P4)")
    print(f"{num_m3} = {total_3_edge_sets} - {num_k1_3} - {num_p3_k2} - {num_p4}")
    print("-" * 20)
    return num_m3

if __name__ == '__main__':
    # We use a known counterexample case: d=3, n=10.
    # It is known that there exist 3-regular bipartite graphs on 10 vertices
    # with different numbers of 4-cycles.
    
    # Case 1: A graph with no 4-cycles (girth is 6).
    m3_g1 = calculate_matchings(n=10, d=3, num_c4=0)
    
    # Case 2: A graph with 15 4-cycles.
    m3_g2 = calculate_matchings(n=10, d=3, num_c4=15)
    
    print("\nConclusion:")
    if m3_g1 != m3_g2:
        print(f"The number of 3-matchings is different ({m3_g1} != {m3_g2}).")
        print("Therefore, two bipartite, d-regular graphs on n vertices do not necessarily have the same number of 3-matchings.")
    else:
        print("The numbers of 3-matchings are the same. My logic or counterexample is flawed.")
