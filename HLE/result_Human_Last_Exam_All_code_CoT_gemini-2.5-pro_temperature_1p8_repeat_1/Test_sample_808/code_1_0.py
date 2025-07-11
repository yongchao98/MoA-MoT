def solve_graph_problem():
    """
    Identifies non-planar graphs after modification from the family of Platonic solids.
    """

    # 1. Define the family of graphs (Platonic solids) with their properties.
    # V is the number of vertices, k is the number of edges for each face.
    platonic_solids = [
        {'name': 'Tetrahedron',   'V': 4,  'k': 3},
        {'name': 'Cube',          'V': 8,  'k': 4},
        {'name': 'Octahedron',    'V': 6,  'k': 3},
        {'name': 'Dodecahedron',  'V': 20, 'k': 5},
        {'name': 'Icosahedron',   'V': 12, 'k': 3}
    ]

    # 2. Order the graphs by their increasing number of vertices.
    sorted_solids = sorted(platonic_solids, key=lambda x: x['V'])

    non_planar_labels = []

    # 3. Iterate through the sorted list, assign labels, and check for planarity.
    for i, solid in enumerate(sorted_solids):
        # Labels are 1-based.
        label = i + 1
        k = solid['k']

        # The modification on a k-gon face creates a K_k subgraph (complete graph on k vertices).
        # A graph containing a K_k subgraph is non-planar if k >= 5.
        # This is the "final equation" or logic to determine the result.
        is_now_non_planar = (k >= 5)

        if is_now_non_planar:
            non_planar_labels.append(str(label))

    # 4. Print the final answer as a comma-delimited list of labels with no spaces.
    final_answer = ",".join(non_planar_labels)
    print(final_answer)

solve_graph_problem()