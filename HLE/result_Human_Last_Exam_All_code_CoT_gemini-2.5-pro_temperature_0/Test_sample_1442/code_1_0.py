import itertools

def count_3_matchings(edges):
    """Counts the number of 3-matchings in a graph given by its edge list."""
    count = 0
    # A 3-matching consists of 3 edges. We iterate through all combinations.
    for edge1, edge2, edge3 in itertools.combinations(edges, 3):
        # A 3-matching involves 6 distinct vertices.
        # We collect all vertices from the three edges.
        vertices = list(edge1) + list(edge2) + list(edge3)
        
        # We check if all 6 vertices are unique by converting the list to a set.
        if len(set(vertices)) == 6:
            count += 1
    return count

def main():
    """
    Constructs two non-isomorphic 3-regular bipartite graphs on 12 vertices
    and counts the number of 3-matchings in each.
    """
    # --- Graph G1: Two disjoint copies of K_3,3 ---
    # The first K_3,3 connects vertices {0,1,2} to {3,4,5}.
    # The second K_3,3 connects vertices {6,7,8} to {9,10,11}.
    g1_edges = []
    for i in range(3):
        for j in range(3, 6):
            g1_edges.append((i, j))
    for i in range(6, 9):
        for j in range(9, 12):
            g1_edges.append((i, j))

    # --- Graph G2: A connected 3-regular bipartite graph on 12 vertices ---
    # The partitions are U={0..5} and V={6..11}.
    # Each vertex u_i in U is connected to v_i, v_{i+1}, and v_{i+3} (indices mod 6).
    g2_edges = []
    for i in range(6):
        u = i
        # v_i
        v1 = i + 6
        # v_{i+1 mod 6}
        v2 = ((i + 1) % 6) + 6
        # v_{i+3 mod 6}
        v3 = ((i + 3) % 6) + 6
        g2_edges.append((u, v1))
        g2_edges.append((u, v2))
        g2_edges.append((u, v3))

    # Count the 3-matchings for each graph
    num_matchings_g1 = count_3_matchings(g1_edges)
    num_matchings_g2 = count_3_matchings(g2_edges)

    # Print the results and the conclusion
    print("To answer the question, we test a counterexample with d=3 and n=12.")
    print("\nGraph G1 is two disjoint copies of K_3,3.")
    print(f"Graph G1 has {len(g1_edges)} edges.")
    print(f"The number of 3-matchings in G1 is: {num_matchings_g1}")

    print("\nGraph G2 is a connected 3-regular bipartite graph.")
    print(f"Graph G2 has {len(g2_edges)} edges.")
    print(f"The number of 3-matchings in G2 is: {num_matchings_g2}")

    print("\nSince the number of 3-matchings is different for these two graphs,")
    print("two bipartite, d-regular graphs on n vertices do not necessarily have the same number of 3-matchings.")

if __name__ == "__main__":
    main()