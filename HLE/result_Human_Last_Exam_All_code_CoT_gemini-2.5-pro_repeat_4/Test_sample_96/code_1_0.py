def solve_e8_torsion_count():
    """
    Calculates the number of minimal length, positive word, order 10 torsion elements in A/Z for type E8.

    This is equivalent to counting the number of parabolic subgroups of type A1 x A4 in the E8 Dynkin diagram.
    """
    # The E8 Dynkin diagram is represented as an adjacency list.
    # The nodes are numbered 1 to 8.
    #     2
    #     |
    # 1 - 3 - 4 - 5 - 6 - 7 - 8
    adj = {
        1: {3},
        2: {3},
        3: {1, 2, 4},
        4: {3, 5},
        5: {4, 6},
        6: {5, 7},
        7: {6, 8},
        8: {7}
    }
    nodes = set(adj.keys())

    # Step 1: Find all subgraphs of type A4 (simple paths of 4 nodes).
    # We use a depth-first search from each node to find paths of length 3 (4 nodes).
    paths = []
    for start_node in nodes:
        stack = [([start_node], {start_node})]  # A tuple of (current_path, visited_nodes)
        while stack:
            path, visited = stack.pop()
            if len(path) == 4:
                # Use a sorted tuple as a canonical representation for the set of nodes in the path
                paths.append(tuple(sorted(path)))
                continue

            last_node = path[-1]
            for neighbor in adj.get(last_node, []):
                if neighbor not in visited:
                    new_path = path + [neighbor]
                    new_visited = visited.copy()
                    new_visited.add(neighbor)
                    stack.append((new_path, new_visited))

    # The set of unique A4 subgraphs (represented by their sorted node tuples)
    unique_a4_subgraphs = sorted(list(set(paths)))

    print("The number of minimal length torsion elements of order 10 is the number of")
    print("parabolic subgroups of type A1 x A4 in the E8 diagram.")
    print("\nWe count these by finding all A4 subgraphs (chains of 4 nodes) and then,")
    print("for each, counting the number of disconnected A1 subgraphs (single nodes).\n")

    total_count = 0
    counts_per_a4 = []

    # Step 2: For each A4 subgraph, count the number of disconnected A1 nodes.
    for a4_nodes_tuple in unique_a4_subgraphs:
        a4_nodes = set(a4_nodes_tuple)

        # Find all neighbors of the A4 subgraph
        neighbors_of_a4 = set()
        for node in a4_nodes:
            neighbors_of_a4.update(adj[node])

        # A1 nodes must not be in the A4 subgraph or be adjacent to it.
        forbidden_nodes = a4_nodes.union(neighbors_of_a4)
        available_a1_nodes = nodes - forbidden_nodes

        count = len(available_a1_nodes)
        total_count += count
        counts_per_a4.append(str(count))
        
        # To make node set printing prettier
        a4_str = str(sorted(list(a4_nodes))).replace('[', '{').replace(']', '}')
        a1_str = str(sorted(list(available_a1_nodes))).replace('[', '{').replace(']', '}')
        print(f"For A4 subgraph on nodes {a4_str}, there are {count} choices for A1: {a1_str}")


    # Step 3: Print the final calculation and result.
    equation_str = " + ".join(counts_per_a4)
    print(f"\nThe total number of such subgroups is the sum of these counts:")
    print(f"{equation_str} = {total_count}")

solve_e8_torsion_count()