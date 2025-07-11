import itertools
from math import comb

def get_edges(adj):
    """Extracts a list of unique edges from an adjacency list."""
    edges = set()
    for u, neighbors in adj.items():
        for v in neighbors:
            if u < v:
                edges.add((u, v))
    return list(edges)

def count_c4_cycles(adj):
    """Counts the number of 4-cycles (C4) in a graph."""
    c4_count = 0
    nodes = list(adj.keys())
    # Iterate over all pairs of distinct vertices
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            u, v = nodes[i], nodes[j]
            # Find common neighbors
            common_neighbors = set(adj[u]) & set(adj[v])
            # Each pair of common neighbors forms a C4 with u and v
            c4_count += comb(len(common_neighbors), 2)
    return c4_count

def brute_force_count_3_matchings(adj):
    """Counts 3-matchings by iterating through all edge combinations."""
    edges = get_edges(adj)
    count = 0
    # Iterate through all combinations of 3 edges
    for edge_combo in itertools.combinations(edges, 3):
        e1, e2, e3 = edge_combo
        # Collect all vertices from the three edges
        vertices = list(e1) + list(e2) + list(e3)
        # A 3-matching has 6 distinct vertices
        if len(set(vertices)) == 6:
            count += 1
    return count

def solve():
    """
    Solves the problem by analyzing two specific graphs G1 and G2,
    calculating their 3-matchings using both a formula and brute force.
    """
    # Graph G1: Tietze's Graph
    # 3-regular, bipartite, 10 vertices
    adj_G1 = {
        0: [1, 5, 9], 1: [0, 2, 6], 2: [1, 3, 7], 3: [2, 4, 8], 4: [3, 5, 9],
        5: [0, 4, 6], 6: [1, 5, 7], 7: [2, 6, 8], 8: [3, 7, 9], 9: [0, 4, 8]
    }

    # Graph G2: Another 3-regular bipartite graph on 10 vertices
    adj_G2 = {
        0: [5, 6, 9], 1: [5, 7, 8], 2: [6, 7, 9], 3: [5, 7, 8], 4: [6, 8, 9],
        5: [0, 1, 3], 6: [0, 2, 4], 7: [1, 2, 3], 8: [1, 3, 4], 9: [0, 2, 4]
    }

    # Graph properties
    n = 10
    d = 3
    num_edges = (n * d) // 2

    # Calculate constant terms (same for both graphs)
    total_triples = comb(num_edges, 3)
    n_k13 = n * comb(d, 3)
    n_p3_k2 = n * comb(d, 2) * (num_edges - (3 * d - 2))

    # --- Calculations for G1 ---
    c4_g1 = count_c4_cycles(adj_G1)
    n_p4_g1 = (num_edges * (d - 1)**2 - 4 * c4_g1) // 2
    m3_g1_formula = total_triples - n_k13 - n_p3_k2 - n_p4_g1
    m3_g1_brute_force = brute_force_count_3_matchings(adj_G1)

    print("--- For Graph G1 ---")
    print(f"Number of 4-cycles c4(G1) = {c4_g1}")
    print("Equation for M3(G1):")
    print(f"M3(G1) = C(|E|,3) - N(K1,3) - N(P3 U K2) - N(P4)")
    print(f"M3(G1) = {total_triples} - {n_k13} - {n_p3_k2} - {n_p4_g1}")
    print(f"M3(G1) = {m3_g1_formula}")
    print(f"Verification by brute force: {m3_g1_brute_force}\n")

    # --- Calculations for G2 ---
    c4_g2 = count_c4_cycles(adj_G2)
    n_p4_g2 = (num_edges * (d - 1)**2 - 4 * c4_g2) // 2
    m3_g2_formula = total_triples - n_k13 - n_p3_k2 - n_p4_g2
    m3_g2_brute_force = brute_force_count_3_matchings(adj_G2)

    print("--- For Graph G2 ---")
    print(f"Number of 4-cycles c4(G2) = {c4_g2}")
    print("Equation for M3(G2):")
    print(f"M3(G2) = C(|E|,3) - N(K1,3) - N(P3 U K2) - N(P4)")
    print(f"M3(G2) = {total_triples} - {n_k13} - {n_p3_k2} - {n_p4_g2}")
    print(f"M3(G2) = {m3_g2_formula}")
    print(f"Verification by brute force: {m3_g2_brute_force}\n")

    if m3_g1_formula != m3_g2_formula:
        print("Conclusion: The number of 3-matchings is NOT necessarily the same.")
    else:
        print("Conclusion: The number of 3-matchings is the same for these two graphs.")

if __name__ == '__main__':
    solve()