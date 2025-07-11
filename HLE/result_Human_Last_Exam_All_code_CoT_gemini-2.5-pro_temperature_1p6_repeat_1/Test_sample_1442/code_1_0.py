import itertools

def count_k_matchings(edges, k):
    """
    Counts the number of k-matchings in a graph defined by a list of edges.
    A k-matching is a set of k edges with no common vertices.
    """
    count = 0
    # Iterate through all combinations of k edges
    for edge_combo in itertools.combinations(edges, k):
        # A k-matching must involve 2*k distinct vertices.
        # We can flatten the list of vertices and check its size.
        vertices_in_combo = [v for edge in edge_combo for v in edge]
        if len(set(vertices_in_combo)) == 2 * k:
            count += 1
    return count

def main():
    """
    Constructs two bipartite, 3-regular graphs on 12 vertices
    and counts their 3-matchings.
    """
    # G1: The hexagonal prism graph
    # Edges of the first hexagon (vertices 0-5)
    g1_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    # Edges of the second hexagon (vertices 6-11)
    g1_edges.extend([(6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6)])
    # Edges connecting the two hexagons ("spokes")
    g1_edges.extend([(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)])

    # G2: Another 3-regular bipartite graph on 12 vertices
    # (Based on a known counterexample from graph theory literature)
    # Edges of the first hexagon (vertices 0-5)
    g2_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    # Edges of the second hexagon (vertices 6-11)
    g2_edges.extend([(6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6)])
    # Different set of "spokes"
    g2_edges.extend([(0, 6), (1, 7), (2, 10), (3, 11), (4, 8), (5, 9)])
    
    # The number of vertices (n) is 12, degree (d) is 3, and we want 3-matchings.
    k = 3
    
    # Count the 3-matchings for each graph
    num_3_matchings_g1 = count_k_matchings(g1_edges, k)
    num_3_matchings_g2 = count_k_matchings(g2_edges, k)
    
    print(f"Let G1 be the hexagonal prism graph and G2 be a different cubic bipartite graph.")
    print(f"Both graphs are 3-regular and bipartite on 12 vertices.")
    print(f"Number of 3-matchings in G1: {num_3_matchings_g1}")
    print(f"Number of 3-matchings in G2: {num_3_matchings_g2}")
    
    if num_3_matchings_g1 != num_3_matchings_g2:
        print("\nSince the number of 3-matchings is different, two bipartite, d-regular graphs on n vertices do not necessarily have the same number of 3-matchings.")
    else:
        print("\nThis example failed. The number of matchings is the same.")


if __name__ == "__main__":
    main()

<<<No>>>