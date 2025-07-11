import numpy as np

def solve():
    """
    Analyzes the given directed graph to determine the correct statement among the choices.
    """
    # Adjacency list for the directed graph (0-indexed)
    # adj[i] = list of vertices that i points to
    adj = {
        0: [3],
        1: [8],
        2: [4],
        3: [2, 6],
        4: [1, 5, 6],
        5: [2, 8],
        6: [1, 5],
        7: [],
        8: [7]
    }
    num_vertices = 9

    # --- Statement A ---
    print("--- Evaluating Statement A ---")
    out_degrees = {v: len(adj[v]) for v in range(num_vertices)}
    max_out_degree_vertex = -1
    max_out_degree_val = -1
    for v, deg in out_degrees.items():
        if deg > max_out_degree_val:
            max_out_degree_val = deg
            max_out_degree_vertex = v
    
    print(f"The vertex with the largest out-degree is {max_out_degree_vertex} with deg+({max_out_degree_vertex}) = {max_out_degree_val}.")
    statement_a_val = 4
    if max_out_degree_vertex == 4 and max_out_degree_val == statement_a_val:
        print("Statement A is correct.")
    else:
        print(f"Statement A claims deg+(4) = {statement_a_val}, which is incorrect.")
    print("-" * 20 + "\n")

    # --- Statement B ---
    print("--- Evaluating Statement B ---")
    walk1 = [0, 3, 2, 4, 6, 1, 8, 7]
    walk2 = [0, 3, 2, 5, 8, 7]

    def is_valid_walk(walk, graph_adj):
        for i in range(len(walk) - 1):
            if walk[i+1] not in graph_adj[walk[i]]:
                return False
        return True

    def calculate_path_length(walk):
        length = 0
        weights = []
        for i in range(len(walk) - 1):
            u, v = walk[i], walk[i+1]
            weight = v - u
            weights.append(weight)
            length += weight
        
        equation_str = " + ".join(f"({w})" for w in weights)
        print(f"Path length = {equation_str} = {length}")
        return length

    print("Checking walk 1: 0->3->2->4->6->1->8->7")
    if is_valid_walk(walk1, adj):
        print("Walk 1 is valid.")
        len1 = calculate_path_length(walk1)
    else:
        print("Walk 1 is invalid.")
        len1 = None

    print("\nChecking walk 2: 0->3->2->5->8->7")
    if is_valid_walk(walk2, adj):
        print("Walk 2 is valid.")
        len2 = calculate_path_length(walk2)
    else:
        print("Walk 2 is invalid because the arc 2->5 does not exist.")
        len2 = None
    
    if len1 is not None and len2 is not None:
        if len1 == len2 + 2:
            print("Statement B is correct.")
        else:
            print("Statement B is incorrect.")
    else:
        print("Statement B is incorrect as at least one walk is invalid.")
    print("-" * 20 + "\n")

    # --- Statement C ---
    print("--- Evaluating Statement C ---")
    given_matrix = np.array([
        [0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0]
    ])
    
    correct_matrix = np.zeros((num_vertices, num_vertices), dtype=int)
    for u, neighbors in adj.items():
        for v in neighbors:
            correct_matrix[u, v] = 1
            
    if np.array_equal(given_matrix, correct_matrix):
        print("Statement C is correct.")
    else:
        print("Statement C is incorrect. The provided matrix does not match the graph.")
        print("Correct adjacency matrix:")
        print(correct_matrix)
        print("Provided matrix:")
        print(given_matrix)
    print("-" * 20 + "\n")

    # --- Statement D ---
    print("--- Evaluating Statement D ---")
    # Part 1: Walks from 3 to 7 in digraph, k <= 10
    M = correct_matrix
    walks_3_7 = 0
    M_k = np.copy(M)
    for k in range(1, 11):
        if k > 1:
            M_k = np.dot(M_k, M)
        walks_3_7 += M_k[3, 7]
    print(f"Number of walks in digraph from 3 to 7 with length k <= 10 is: {walks_3_7}")

    # Part 2: Walks from 6 to 5 in undirected graph, k <= 3
    U = np.zeros((num_vertices, num_vertices), dtype=int)
    for i in range(num_vertices):
        for j in range(num_vertices):
            if M[i, j] == 1 or M[j, i] == 1:
                U[i, j] = 1
    
    walks_6_5 = 0
    U_k = np.copy(U)
    for k in range(1, 4):
        if k > 1:
            U_k = np.dot(U_k, U)
        walks_6_5 += U_k[6, 5]
    print(f"Number of walks in undirected graph from 6 to 5 with length k <= 3 is: {walks_6_5}")

    if walks_3_7 == walks_6_5:
        print(f"The numbers are equal ({walks_3_7} == {walks_6_5}). Statement D is correct.")
    else:
        print(f"The numbers are not equal ({walks_3_7} != {walks_6_5}). Statement D is incorrect.")
    print("-" * 20 + "\n")

    # --- Statement E ---
    print("--- Evaluating Statement E ---")
    # The "degree sum" in this context likely refers to the sum of out-degrees (or in-degrees), which equals the number of arcs.
    num_arcs = sum(len(v) for v in adj.values())
    degree_list = [len(adj[v]) for v in range(num_vertices)]
    equation_str = " + ".join(map(str, degree_list))
    print(f"The sum of out-degrees is the total number of arcs: {equation_str} = {num_arcs}")
    if num_arcs == 13:
        print("Statement E is correct.")
    else:
        print("Statement E is incorrect.")
    print("-" * 20 + "\n")

solve()