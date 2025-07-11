import itertools

def has_c4(adj_list):
    """
    Checks if a graph represented by an adjacency list contains a C4 (quadrilateral).
    A C4 is a cycle of length 4: a-b-c-d-a.
    This can be found by looking for two vertices (a, c) that share two distinct common neighbors (b, d).
    """
    nodes = list(adj_list.keys())
    # Iterate over all distinct pairs of nodes (a, c)
    for a, c in itertools.combinations(nodes, 2):
        # Find common neighbors of a and c
        # Use sets for efficient intersection
        neighbors_of_a = set(adj_list.get(a, []))
        neighbors_of_c = set(adj_list.get(c, []))
        
        common_neighbors = list(neighbors_of_a.intersection(neighbors_of_c))
        
        # If there are 2 or more common neighbors, a C4 exists.
        # For example, if b and d are common neighbors, then a-b-c-d-a is a C4.
        if len(common_neighbors) >= 2:
            b, d = common_neighbors[0], common_neighbors[1]
            print(f"  - Found C4: {a}-{b}-{c}-{d}-{a}")
            return True
    return False

def main():
    """
    Main function to construct the graph, check for C4, and test for maximality.
    """
    # Adjacency list for the constructed graph on 8 vertices {0, 1, ..., 7}
    # Vertex 0 is the central vertex 'C'
    # Vertices 1-7 are the outer vertices 'v1' to 'v7'
    graph = {
        0: [1, 2, 3, 4, 5, 6, 7],  # Central vertex connected to all others
        1: [0, 2],
        2: [0, 1],
        3: [0, 4],
        4: [0, 3],
        5: [0, 6],
        6: [0, 5],
        7: [0]
    }
    
    num_edges = sum(len(v) for v in graph.values()) // 2
    print(f"Constructed graph has {len(graph)} vertices and {num_edges} edges.")

    print("\nStep 1: Checking if the constructed graph has a C4...")
    if not has_c4(graph):
        print("  - The graph is C4-free.")
    else:
        print("  - The graph contains a C4.")

    print("\nStep 2: Checking if the graph is a maximal C4-free graph.")
    print("This is done by trying to add every possible non-edge and checking if it creates a C4.")
    
    is_maximal = True
    nodes = list(graph.keys())
    
    # Generate all possible pairs of vertices (potential edges)
    for u, v in itertools.combinations(nodes, 2):
        # Check if the edge (u, v) already exists
        if v not in graph[u]:
            # This is a non-edge. Let's try to add it.
            
            # Create a copy of the graph to avoid modifying the original
            temp_graph = {node: neighbors[:] for node, neighbors in graph.items()}
            
            # Add the new edge
            temp_graph[u].append(v)
            temp_graph[v].append(u)
            
            # Check if this new edge creates a C4
            print(f"Testing addition of edge ({u}, {v})...")
            if not has_c4(temp_graph):
                is_maximal = False
                print(f"  - Adding edge ({u}, {v}) did NOT create a C4.")
                # We can stop as soon as we find one, but let's check all for completeness.

    if is_maximal:
        print("\nConclusion: The constructed graph is maximal. Adding any single edge creates a C4.")
    else:
        print("\nConclusion: The constructed graph is not maximal.")

    # Final Answer
    final_answer = num_edges
    print(f"\nThe equation representing the number of edges is: 7 (star edges) + 3 (pairing edges) = 10")
    print(f"The maximum number of edges in a C4-free graph with 8 vertices is {final_answer}.")


if __name__ == "__main__":
    main()
