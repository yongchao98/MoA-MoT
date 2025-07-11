import itertools

def count_k_matchings(graph, k):
    """
    Counts the number of k-matchings in a graph using a backtracking algorithm.
    
    A k-matching is a set of k edges where no two edges share a vertex.
    
    Args:
        graph (dict): An adjacency list representation of the graph.
        k (int): The size of the matching.
        
    Returns:
        int: The number of k-matchings.
    """
    edges = set()
    for u, neighbors in graph.items():
        for v in neighbors:
            # Store edges uniquely, e.g., (min(u,v), max(u,v))
            if u < v:
                edges.add((u, v))
    
    edges = list(edges)
    n_edges = len(edges)
    count = 0
    
    # Use itertools.combinations to find all sets of k edges
    for combo in itertools.combinations(edges, k):
        vertices_used = set()
        is_matching = True
        for edge in combo:
            if edge[0] in vertices_used or edge[1] in vertices_used:
                is_matching = False
                break
            vertices_used.add(edge[0])
            vertices_used.add(edge[1])
        
        if is_matching:
            count += 1
            
    return count

def define_graphs():
    """
    Defines the two graphs G1 (Y_6) and G2 (K_3,3 U K_3,3)
    """
    # G1: Prism graph Y_6 = C_6 x K_2
    # Vertices 0-5 form one C_6, vertices 6-11 form the other.
    # Edges are (i, i+1) on each cycle and spokes (i, i+6)
    G1 = {i: set() for i in range(12)}
    for i in range(6):
        # Rim edges
        G1[i].add((i + 1) % 6)
        G1[(i + 1) % 6].add(i)
        G1[i + 6].add((i + 1) % 6 + 6)
        G1[(i + 1) % 6 + 6].add(i + 6)
        # Spoke edges
        G1[i].add(i + 6)
        G1[i + 6].add(i)
        
    # G2: Disjoint union of two K_3,3 graphs
    # First K_3,3 on vertices 0-5. Partitions U={0,1,2}, W={3,4,5}
    # Second K_3,3 on vertices 6-11. Partitions U={6,7,8}, W={9,10,11}
    G2 = {i: set() for i in range(12)}
    # First component
    for i in range(3):
        for j in range(3, 6):
            G2[i].add(j)
            G2[j].add(i)
    # Second component
    for i in range(6, 9):
        for j in range(9, 12):
            G2[i].add(j)
            G2[j].add(i)

    return G1, G2

def main():
    G1, G2 = define_graphs()
    k = 3
    
    num_matchings_g1 = count_k_matchings(G1, k)
    num_matchings_g2 = count_k_matchings(G2, k)
    
    print(f"Let G1 be the prism graph Y_6 and G2 be the disjoint union of two K_3,3 graphs.")
    print(f"Both graphs are 3-regular, bipartite, and have 12 vertices.")
    print(f"Number of 3-matchings in G1 (Y_6): {num_matchings_g1}")
    print(f"Number of 3-matchings in G2 (K_3,3 U K_3,3): {num_matchings_g2}")
    
    if num_matchings_g1 != num_matchings_g2:
        print("\nThe number of 3-matchings is different.")
        print("Therefore, two bipartite, d-regular graphs on n vertices do not necessarily have the same number of 3-matchings.")
    else:
        print("\nThe number of 3-matchings is the same for this pair of graphs.")

if __name__ == "__main__":
    main()
