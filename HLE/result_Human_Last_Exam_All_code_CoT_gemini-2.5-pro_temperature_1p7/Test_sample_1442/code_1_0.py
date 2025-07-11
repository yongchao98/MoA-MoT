import itertools

def count_3_matchings(adj):
    """Counts the number of 3-matchings in a graph given its adjacency list."""
    
    # Step 1: Generate the list of edges from the adjacency list.
    # We use a set to store edges to avoid duplicates like (u,v) and (v,u).
    # Each edge is stored as a sorted tuple to ensure uniqueness.
    edges = set()
    for u, neighbors in adj.items():
        for v in neighbors:
            edges.add(tuple(sorted((u, v))))
    
    edge_list = list(edges)
    num_edges = len(edge_list)
    
    # Check if there are enough edges to form a 3-matching
    if num_edges < 3:
        return 0

    # Step 2: Iterate through all combinations of 3 edges.
    count = 0
    for edge_combo in itertools.combinations(edge_list, 3):
        # edge_combo is a tuple of 3 edges, e.g., ((u1, v1), (u2, v2), (u3, v3))
        
        # Step 3: Check for disjointness.
        # Flatten the list of vertices from the 3 edges
        vertices_in_combo = [v for edge in edge_combo for v in edge]
        
        # A 3-matching has 6 unique vertices.
        if len(set(vertices_in_combo)) == 6:
            count += 1
            
    return count

def main():
    # Graph G1: The Cube Graph (Q_3)
    # Vertices are mapped from binary strings {000, ..., 111} to integers {0, ..., 7}.
    # Edge exists if Hamming distance is 1.
    g1_adj = {
        0: [1, 2, 4], 1: [0, 3, 5], 2: [0, 3, 6], 3: [1, 2, 7],
        4: [0, 5, 6], 5: [1, 4, 7], 6: [2, 4, 7], 7: [3, 5, 6]
    }

    # Graph G2: The non-cubical bipartite 3-regular graph on 8 vertices.
    # It's the wreath product of C4 and K2.
    g2_adj = {
        0: [4, 5, 6], 1: [5, 6, 7], 2: [4, 6, 7], 3: [4, 5, 7],
        4: [0, 2, 3], 5: [0, 1, 3], 6: [0, 1, 2], 7: [1, 2, 3]
    }

    # Both graphs are bipartite, 3-regular, and have 8 vertices.
    # Total number of edges in a 3-regular graph on 8 vertices is (8 * 3) / 2 = 12.
    
    # Count 3-matchings for each graph
    g1_3_matchings = count_3_matchings(g1_adj)
    g2_3_matchings = count_3_matchings(g2_adj)

    print(f"Graph G1 is the Cube graph Q_3.")
    print(f"Number of 3-matchings in G1: {g1_3_matchings}")
    print(f"\nGraph G2 is the other 3-regular bipartite graph on 8 vertices (C4 wreath K2).")
    print(f"Number of 3-matchings in G2: {g2_3_matchings}")

    if g1_3_matchings == g2_3_matchings:
        print("\nThe number of 3-matchings is the same for these two graphs.")
        print("Based on this example, one cannot conclude that the answer is No.")
    else:
        print("\nThe number of 3-matchings is different for these two graphs.")
        print("This serves as a counterexample.")

if __name__ == '__main__':
    main()
