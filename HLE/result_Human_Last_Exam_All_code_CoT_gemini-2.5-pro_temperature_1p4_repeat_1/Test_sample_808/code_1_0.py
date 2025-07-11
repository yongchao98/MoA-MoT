def find_non_planar_graphs():
    """
    This function identifies a class of graphs, orders them, and determines
    which become non-planar after a specified modification.

    The graphs are the five Platonic solids. The modification turns a face
    of degree 'f' into a complete graph K_f. The resulting graph is non-planar
    if f >= 5, as it will contain a K_5 subgraph.
    """

    # Step 1 & 2: Define the Platonic solids with their number of vertices (V)
    # and face degree (f), then order them by V.
    # Format: (Name, V, f)
    platonic_solids = [
        ('Tetrahedron', 4, 3),
        ('Octahedron', 6, 3),
        ('Cube', 8, 4),
        ('Icosahedron', 12, 3),
        ('Dodecahedron', 20, 5)
    ]

    # The list is already sorted by V. We can assign labels.
    # Label is the 1-based index.
    
    # Step 3, 4, 5: Analyze each graph and collect labels of non-planar results.
    non_planar_labels = []
    for i, (name, V, f) in enumerate(platonic_solids):
        label = i + 1
        
        # A face of degree f is turned into a K_f.
        # If f >= 5, the graph contains a K_5 subgraph and is non-planar.
        if f >= 5:
            non_planar_labels.append(str(label))

    # Print the final result as a comma-separated string.
    print(",".join(non_planar_labels))

find_non_planar_graphs()