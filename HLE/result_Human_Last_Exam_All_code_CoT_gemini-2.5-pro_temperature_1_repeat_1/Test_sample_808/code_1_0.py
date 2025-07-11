def find_non_planar_graphs():
    """
    This script identifies which Platonic solids become non-planar after a specific modification.

    The plan is as follows:
    1. Define the five Platonic solids with their number of vertices (V) and face degree (s).
    2. Order them by V to establish their labels (1 to 5).
    3. Analyze the modification: adding edges to a face's boundary to form a complete graph K_s.
    4. A graph becomes non-planar if this procedure creates a K_s subgraph where s >= 5.
    5. Collect the labels of the graphs that meet this criterion.
    """

    # Step 1: Define the properties of the Platonic solids.
    # name: Name of the solid
    # V: Number of vertices
    # s: Number of edges/vertices per face
    platonic_solids = [
        {'name': 'Tetrahedron', 'V': 4, 's': 3},
        {'name': 'Octahedron', 'V': 6, 's': 3},
        {'name': 'Cube', 'V': 8, 's': 4},
        {'name': 'Icosahedron', 'V': 12, 's': 3},
        {'name': 'Dodecahedron', 'V': 20, 's': 5},
    ]

    # Step 2: Order the graphs by increasing number of vertices to assign labels.
    sorted_solids = sorted(platonic_solids, key=lambda x: x['V'])

    non_planar_labels = []

    # Step 3 & 4: Iterate through the sorted graphs and check the planarity condition.
    for i, graph_info in enumerate(sorted_solids):
        label = i + 1
        s = graph_info['s']
        
        # If a face has 5 or more vertices (s >= 5), turning it into a complete
        # graph (K_s) will introduce a K_5 subgraph, making the entire graph non-planar.
        # For s < 5, the resulting graph remains planar.
        if s >= 5:
            non_planar_labels.append(str(label))

    # Step 5: Print the final result in the required format.
    result = ",".join(non_planar_labels)
    print(result)

find_non_planar_graphs()