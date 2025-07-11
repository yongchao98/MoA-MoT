def solve_graph_problem():
    """
    This function analyzes the given directed graph to verify statement E.
    It calculates the sum of the out-degrees of all vertices, which is
    equivalent to the total number of edges in the graph.
    """
    # Adjacency list representation of the directed graph
    # graph[i] contains a list of vertices that vertex i has an edge to.
    graph = {
        0: [3],
        1: [8],
        2: [4],
        3: [2, 6],
        4: [1, 2, 5],
        5: [8],
        6: [1, 5],
        7: [3],
        8: [7]
    }

    out_degrees = {}
    total_degree_sum = 0
    
    print("Calculating the sum of out-degrees (which equals the number of edges):")
    
    # Sort vertices for consistent output order
    sorted_vertices = sorted(graph.keys())
    
    degree_values = []
    for v in sorted_vertices:
        degree = len(graph.get(v, []))
        out_degrees[v] = degree
        total_degree_sum += degree
        degree_values.append(str(degree))
        print(f"deg+({v}) = {degree}")

    equation_str = " + ".join(degree_values)
    print(f"\nTotal degree sum = {equation_str} = {total_degree_sum}")

solve_graph_problem()
<<<E>>>