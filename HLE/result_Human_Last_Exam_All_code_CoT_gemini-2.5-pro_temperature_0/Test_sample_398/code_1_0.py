import networkx as nx

def demonstrate_non_unique_crossing():
    """
    This function demonstrates that adding an edge to a maximal planar graph
    can result in a graph with a non-unique 1-crossing drawing.
    """
    # 1. Define a maximal planar graph G on 5 vertices.
    # This graph is a wheel graph (vertex 0 connected to 1,2,3,4)
    # with the outer cycle (1-2-3-4-1) triangulated by adding the edge (1,3).
    G = nx.Graph()
    vertices = [0, 1, 2, 3, 4]
    G.add_nodes_from(vertices)
    
    # Edges of the maximal planar graph G
    # It has 5 vertices and 3*5-6 = 9 edges.
    edges_G = [
        (0, 1), (0, 2), (0, 3), (0, 4),  # Wheel spokes
        (1, 2), (2, 3), (3, 4), (4, 1),  # Wheel rim
        (1, 3)                           # Triangulating chord
    ]
    G.add_edges_from(edges_G)

    # 2. Define an edge 'e' that is not in G.
    # The only missing edge besides e is (1,3) which is already in G.
    # The only non-adjacent pair is (2, 4).
    e = (2, 4)
    
    print(f"Consider the maximal planar graph G with vertices {G.nodes()} and edges {G.edges()}.")
    print(f"We add the edge e = {e}, which is not in G.")
    print("A 1-crossing drawing of G + e can be made by having edge 'e' cross an existing edge 'f' from G.")
    print("This is possible if the graph G - {f} + {e} is planar.")
    print("\nLet's find all edges 'f' in G that 'e' can cross:")

    possible_crossings = []
    # 3. Iterate through all edges f in G.
    for f in G.edges():
        # 4. Create the graph H = G - {f} + {e}
        H = G.copy()
        H.remove_edge(*f)
        H.add_edge(*e)
        
        # 5. Check if H is planar.
        # networkx's check_planarity returns a tuple (is_planar, certificate)
        # The certificate is a Kuratowski subgraph if not planar.
        is_planar, _ = nx.check_planarity(H)
        
        if is_planar:
            possible_crossings.append(f)

    # 6. Print the results.
    if len(possible_crossings) > 1:
        print(f"\nFound {len(possible_crossings)} different edges in G that the new edge {e} can cross:")
        for f in possible_crossings:
            print(f" - Edge {f}")
        print("\nSince there is more than one choice for the crossed edge, the 1-crossing drawing of G' is not unique.")
    elif len(possible_crossings) == 1:
        print(f"\nFound only one edge {possible_crossings[0]} that can be crossed. The drawing is unique in this case.")
    else:
        print("\nCould not find any edge f to cross to create a 1-crossing drawing.")

demonstrate_non_unique_crossing()