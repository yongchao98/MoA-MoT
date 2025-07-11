def solve():
    """
    This script verifies Statement A by analyzing the given directed graph.
    """
    # The directed graph is represented as an adjacency list.
    graph = {
        0: [3],
        1: [8],
        2: [4],
        3: [2, 6],
        4: [1, 2, 5, 6],
        5: [8],
        6: [1, 5],
        7: [],
        8: [7]
    }
    
    print("Evaluating Statement A: The vertex v with the largest value of deg^+(v) is the vertex labeled 4, with deg^+(4)=4.\n")

    # Step 1: Calculate the out-degree for vertex 4.
    # The out-degree deg^+(v) is the number of outgoing edges from vertex v.
    # The "equation" for the out-degree of vertex 4 is deg^+(4) = |{neighbors of 4}|.
    v = 4
    neighbors_of_v = graph[v]
    degree_of_v = len(neighbors_of_v)
    
    print(f"For vertex {v}, the outgoing edges go to vertices: {neighbors_of_v}")
    print(f"The out-degree is the count of these neighbors.")
    print(f"deg^+({v}) = |{neighbors_of_v}| = {degree_of_v}")

    # Step 2: Find the maximum out-degree among all vertices to confirm it's the largest.
    out_degrees = {vertex: len(neighbors) for vertex, neighbors in graph.items()}
    max_out_degree_vertex = -1
    max_out_degree_value = -1
    
    for vertex, degree in out_degrees.items():
        if degree > max_out_degree_value:
            max_out_degree_value = degree
            max_out_degree_vertex = vertex

    print(f"\nThe maximum out-degree found in the graph is {max_out_degree_value}, for vertex {max_out_degree_vertex}.")

    # Step 3: Compare with the statement.
    print("\n--- Conclusion ---")
    if max_out_degree_vertex == 4 and max_out_degree_value == 4:
        print("The calculation matches the statement.")
        print("Statement A is correct.")
    else:
        print("The calculation does not match the statement.")
        print("Statement A is incorrect.")

solve()