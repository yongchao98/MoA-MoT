import collections

def solve_graph_problem():
    """
    This function identifies the graphs in question, analyzes the specified
    modification, and determines which of the resulting graphs are non-planar.
    """
    
    # Step 1: Define the family of graphs.
    # The conditions describe the graphs of the five Platonic solids.
    # We define them by name, number of vertices (V), and vertices per face (k).
    Graph = collections.namedtuple('Graph', ['name', 'V', 'k'])
    platonic_solids = [
        Graph(name='Tetrahedron', V=4, k=3),
        Graph(name='Cube', V=8, k=4),
        Graph(name='Octahedron', V=6, k=3),
        Graph(name='Dodecahedron', V=20, k=5),
        Graph(name='Icosahedron', V=12, k=3),
    ]

    # Step 2: Order the graphs by vertex count and assign labels.
    ordered_graphs = sorted(platonic_solids, key=lambda g: g.V)
    
    # Add the label to each graph object
    LabelledGraph = collections.namedtuple('LabelledGraph', ['label', 'name', 'V', 'k'])
    final_graph_list = []
    for i, g in enumerate(ordered_graphs):
        final_graph_list.append(LabelledGraph(label=i + 1, name=g.name, V=g.V, k=g.k))
        
    print("Analyzing the specified family of graphs:")

    non_planar_labels = []

    # Step 3 & 4: Analyze each graph.
    # The modification turns a face with k vertices into a K_k subgraph.
    # The resulting graph is non-planar if this subgraph is K_5 or larger.
    # The key condition is whether k >= 5.
    
    for graph in final_graph_list:
        print(f"\n--- Processing Graph {graph.label}: The {graph.name} ---")
        print(f"This graph has regions with k = {graph.k} vertices.")

        if graph.k >= 5:
            # This is the condition for non-planarity.
            # Here we "output each number in the final equation!" by explaining the k=5 check.
            print(f"The procedure transforms one face into a K_{graph.k} subgraph.")
            print(f"Because a graph containing a K_5 subgraph is non-planar, and k={graph.k} is >= 5, the modified graph is non-planar.")
            non_planar_labels.append(str(graph.label))
        elif graph.k == 4:
            # Case for the Cube
            print(f"The procedure transforms one face into a K_4 subgraph.")
            print(f"This makes one square face of the cube into a K_4, but the overall graph remains planar.")
        else: # graph.k == 3
            # Case for triangles
            print(f"The regions are triangles (k=3), which have no non-adjacent vertices.")
            print(f"No edges are added, so the graph remains planar.")
    
    # Step 5: Consolidate and print the final answer.
    result = ",".join(sorted(non_planar_labels, key=int))
    
    print("\n-----------------------------------------------------")
    print("The final list of labels for the non-planar graphs is:")
    print(result)

# Execute the solution
solve_graph_problem()