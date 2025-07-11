import collections

def reduce_disjoint_triangles_to_disjoint_cycles(g_tri, k_tri):
    """
    Performs a reduction from k-Disjoint Triangles to DisjointCycles.

    Args:
        g_tri (dict): The adjacency list of the graph for k-Disjoint Triangles.
        k_tri (int): The number of disjoint triangles to find.
    
    Prints:
        The adjacency list of the new graph for the DisjointCycles problem.
        The new parameter k for the DisjointCycles problem.
    """
    if k_tri < 1:
        print("k must be at least 1.")
        return

    # The parameter for DisjointCycles is the same as for k-Disjoint Triangles.
    k_dc = k_tri
    
    # Number of new vertices to add to each edge.
    # To make a triangle (length 3) into a cycle of length at least k_dc,
    # we need to lengthen each of its 3 edges.
    # If we subdivide each edge s times, the new length is 3 + 3s.
    # We need 3 + 3s >= k_dc. Let's choose s = k_dc - 1.
    # The new length will be 3 + 3(k_dc - 1) = 3k_dc, which is >= k_dc.
    subdivisions = k_dc - 1

    g_dc = collections.defaultdict(list)
    
    # Keep track of processed edges to avoid duplicating subdivisions
    processed_edges = set()

    # Add original vertices to the new graph
    for u in g_tri:
        g_dc[u] # ensures vertex exists in the dict

    # Iterate over all edges and subdivide them
    for u in g_tri:
        for v in g_tri[u]:
            # Process each edge only once (e.g., from smaller to larger vertex)
            edge = tuple(sorted((u, v)))
            if edge in processed_edges:
                continue
            
            processed_edges.add(edge)
            
            last_vertex = u
            
            # Create a path of new vertices between u and v
            for i in range(subdivisions):
                # New vertices are named as tuples to avoid collision with original names
                new_vertex = (u, v, i) if u < v else (v, u, i)
                g_dc[last_vertex].append(new_vertex)
                g_dc[new_vertex].append(last_vertex)
                last_vertex = new_vertex
            
            # Connect the end of the path to v
            g_dc[last_vertex].append(v)
            g_dc[v].append(last_vertex)

    print("--- Reduction Output ---")
    print(f"Original k_tri = {k_tri}")
    print(f"New k_dc for DisjointCycles = {k_dc}")
    print("\nAdjacency list of the new graph G':")
    for node, neighbors in sorted(g_dc.items()):
        # Sort neighbors for consistent output
        sorted_neighbors = sorted(neighbors)
        print(f"  {node}: {sorted_neighbors}")
    print("------------------------")


# Example Usage:
# A graph with two disjoint triangles: (0,1,2) and (3,4,5)
# Let's find 2 disjoint triangles (k_tri = 2)
g_example = {
    0: [1, 2],
    1: [0, 2],
    2: [0, 1],
    3: [4, 5],
    4: [3, 5],
    5: [3, 4]
}
k_example = 2
reduce_disjoint_triangles_to_disjoint_cycles(g_example, k_example)
