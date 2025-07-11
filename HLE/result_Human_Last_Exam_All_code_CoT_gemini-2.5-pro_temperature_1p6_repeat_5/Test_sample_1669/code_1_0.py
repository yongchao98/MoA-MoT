import collections

def demonstrate_k_vector():
    """
    This function demonstrates a valid k-vector for a sample graph.

    The main problem asks for the smallest k for ANY 20-vertex, bridgeless,
    3-regular graph. The answer to that question is k=4.

    This is because:
    1. A valid k=3 vector exists if and only if the graph is 3-edge-colorable.
    2. There exist 20-vertex bridgeless 3-regular graphs (called snarks) that are
       NOT 3-edge-colorable. For these graphs, k must be greater than 3.
    3. It has been proven that a valid k=4 vector exists for ALL bridgeless
       cubic graphs.

    Therefore, k=4 is the smallest integer that works for all cases.

    This code will demonstrate the principle by constructing a valid k=3 vector
    for a simpler, 3-edge-colorable graph (a triangular prism) to show what
    such a vector looks like.
    """

    # Define a simple 3-regular graph: the triangular prism.
    # Vertices are 0-5.
    # Edges are tuples of vertices. We use frozenset for undirected edges.
    edges = {
        frozenset({0, 1}), frozenset({1, 2}), frozenset({2, 0}),  # Top face
        frozenset({3, 4}), frozenset({4, 5}), frozenset({5, 3}),  # Bottom face
        frozenset({0, 3}), frozenset({1, 4}), frozenset({2, 5})   # Side edges
    }
    num_vertices = 6

    # For a 3-edge-colorable graph, we can construct a valid 3-vector.
    # We assign values based on a 3-edge-coloring.
    # Let color 1 -> 1, color 2 -> 1, color 3 -> -2.
    # The values {1, -2} are in the set for k=3: {-2, -1, 1, 2}.
    k = 3
    
    # 1. Define a 3-edge-coloring (3 perfect matchings)
    coloring = {
        # Color "Red"
        frozenset({0, 1}), frozenset({2, 5}), frozenset({3, 4}),
        # Color "Green"
        frozenset({1, 2}), frozenset({0, 3}), frozenset({4, 5}),
        # Color "Blue"
        frozenset({2, 0}), frozenset({1, 4}), frozenset({5, 3})
    }
    
    # 2. Assign values to edges based on color to create the k-vector.
    k_vector = {}
    color_to_value = {"Red": 1, "Green": 1, "Blue": -2}
    for color, edge_set in coloring.items():
        for edge in edge_set:
            k_vector[edge] = color_to_value[color]

    print(f"Demonstrating a valid k-vector for k={k} on a triangular prism graph.")
    print(f"Edge values are assigned from the set {{+/-1, ..., +/-{k-1}}}}.")
    print("Value assignment based on a 3-edge-coloring: Red edges=1, Green edges=1, Blue edges=-2\n")

    # 3. Verify the null space condition at each vertex.
    # For each vertex, find its incident edges and sum their values.
    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(frozenset({u, v}))
        adj[v].append(frozenset({u, v}))
        
    print("Verifying the sum of edge values at each vertex is 0:")
    for v in range(num_vertices):
        incident_edges = adj[v]
        vals = [k_vector[edge] for edge in incident_edges]
        edge_names = [str(tuple(e)) for e in incident_edges]
        total = sum(vals)
        
        # Output the equation for each vertex
        # Example: For vertex 0: val(0,1) + val(0,2) + val(0,3) = 1 + (-2) + 1 = 0
        val_strings = [f"{v: >2}" for v in vals]
        expr_vals = " + ".join(val_strings).replace("+ -", "- ")
        edge_str = ' + '.join([f"x{name}" for name in edge_names])
        print(f"Vertex {v}: {edge_str} = {expr_vals} = {total}")

    print("\nAs shown, k=3 works for this graph.")
    print("However, for a non-3-edge-colorable graph (a snark), k=3 would fail.")
    print("The smallest value that works for ALL such graphs is k=4.")


demonstrate_k_vector()