def analyze_decision_coverage():
    """
    Analyzes the decision coverage of a control flow graph based on given test cases.
    """
    # Step 1 & 2: Define decision points and their branches from the graph.
    decision_points = {
        'A': {'A->B', 'A->C'},
        'D': {'D->E', 'D->F'}
    }
    print("Decision Points and Their Branches:")
    for node, branches in decision_points.items():
        print(f"Node {node}: {', '.join(branches)}")
    print("-" * 40)

    # Step 3: Define the executed test paths and extract the edges (branches) they cover.
    test_paths = [
        ['A', 'B', 'D', 'E', 'G'],
        ['A', 'B', 'D', 'E', 'F', 'G'],
        ['A', 'C', 'F', 'C', 'F', 'C', 'F', 'G']
    ]

    covered_branches = set()
    for i, path in enumerate(test_paths, 1):
        path_edges = set()
        for j in range(len(path) - 1):
            edge = f"{path[j]}->{path[j+1]}"
            path_edges.add(edge)
        print(f"Test_{i} path: {', '.join(path)}")
        print(f"Branches covered by Test_{i}: {', '.join(sorted(list(path_edges))) if path_edges else 'None'}")
        covered_branches.update(path_edges)
    print("-" * 40)

    # Step 4: Calculate and report the final coverage.
    print("Coverage Analysis Results:")
    final_conclusion = ""
    all_covered = True

    for node, branches in decision_points.items():
        uncovered = branches - covered_branches
        if not uncovered:
            print(f"[COVERED]   Decision at node {node} is fully covered.")
        else:
            all_covered = False
            print(f"[UNCOVERED] Decision at node {node} is NOT fully covered.")
            print(f"            Missing branches: {', '.join(sorted(list(uncovered)))}")

    print("-" * 40)

    # Step 5: Evaluate the options.
    if all_covered:
        final_conclusion = "B. 100% decision coverage has been achieved."
    else:
        # Check why coverage failed
        if 'D->F' in str(decision_points['D'] - covered_branches):
            final_conclusion = "A. The decision at node D has not been fully tested."
        elif 'E' in decision_points and decision_points['E'] - covered_branches:
             final_conclusion = "C. The decision at node E has not been fully tested."
        elif 'F' in decision_points and decision_points['F'] - covered_branches:
             final_conclusion = "D. The decision at node F has not been fully tested."

    print(f"Final Conclusion: {final_conclusion}")


if __name__ == "__main__":
    analyze_decision_coverage()