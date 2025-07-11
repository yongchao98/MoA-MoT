import itertools

def solve():
    """
    Analyzes a specific maximal planar graph to determine properties of adding a new edge.
    """
    # Define a non-symmetric maximal planar graph G on 5 vertices.
    # This graph is constructed by taking a K4 on {1,2,3,4} and placing vertex 5
    # inside the face (1,2,4) and connecting it to 1, 2, and 4.
    edges = {
        (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4),  # K4 edges
        (1, 5), (2, 5), (4, 5)  # Edges connecting vertex 5
    }
    num_vertices = 5
    num_edges = len(edges)
    
    # A graph is maximal planar if E = 3V - 6
    is_maximal_planar = (num_edges == 3 * num_vertices - 6)
    print(f"The graph G has {num_vertices} vertices and {num_edges} edges.")
    print(f"For a maximal planar graph, expected edges = 3*V - 6 = {3*num_vertices - 6}.")
    if not is_maximal_planar:
        print("The defined graph is not maximal planar.")
        return
        
    print("The graph G is maximal planar.")
    
    # Find a pair of non-adjacent vertices (u, v) to form the new edge e.
    # By construction, (3, 5) is a non-edge.
    u, v = 3, 5
    print(f"\nLet's add the edge e = ({u}, {v}), which is not in G.")

    # Find common neighbors of u and v
    adj = {i: set() for i in range(1, num_vertices + 1)}
    for p1, p2 in edges:
        adj[p1].add(p2)
        adj[p2].add(p1)

    neighbors_u = adj[u]
    neighbors_v = adj[v]
    common_neighbors = neighbors_u.intersection(neighbors_v)
    
    print(f"Neighbors of {u}: {sorted(list(neighbors_u))}")
    print(f"Neighbors of {v}: {sorted(list(neighbors_v))}")
    print(f"Common neighbors of {u} and {v}: {sorted(list(common_neighbors))}")

    # A 1-crossing drawing can be made by crossing an edge f=(x,y) where
    # x and y are common neighbors of u and v.
    # Find all such possible edges f.
    possible_crossed_edges = []
    for x, y in itertools.combinations(common_neighbors, 2):
        # Ensure the pair is sorted to match the set representation
        if tuple(sorted((x,y))) in edges or tuple(sorted((y,x))) in edges:
            possible_crossed_edges.append(tuple(sorted((x,y))))

    print(f"\nTo create a 1-crossing drawing of G + e, the edge e=({u},{v}) can cross any of the following edges f:")
    for edge in possible_crossed_edges:
        print(f" - {edge}")
        
    print("\nSince there are multiple distinct choices for the edge f to be crossed,")
    print("and the graph G is not perfectly symmetric with respect to these choices,")
    print("the resulting 1-crossing drawing is not unique.")
    print("\nThis supports the conclusion that G' can be drawn with at most one crossing, but this drawing is not unique.")

solve()