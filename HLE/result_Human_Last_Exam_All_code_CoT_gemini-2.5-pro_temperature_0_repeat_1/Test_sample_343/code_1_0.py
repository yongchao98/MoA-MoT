def add_universal_vertex(graph):
    """
    Takes a graph represented as a dictionary (adjacency list)
    and returns a new graph with an added universal vertex.
    The new vertex is named 'u*' to avoid collisions.
    """
    new_graph = {v: set(neighbors) for v, neighbors in graph.items()}
    new_vertex = 'u*'
    while new_vertex in new_graph:
        new_vertex += '*'
    
    existing_vertices = set(new_graph.keys())
    new_graph[new_vertex] = existing_vertices
    for vertex in existing_vertices:
        new_graph[vertex].add(new_vertex)
            
    return new_graph

def count_k_cliques_using_oracle(graph, k, p_count_oracle):
    """
    Computes the number of k-cliques in a graph using a hypothetical
    FPT oracle for the PCount problem. This demonstrates the FPT-Turing 
    reduction from #k-Clique to PCount.
    """
    if k < 1:
        return 0
    if k == 1:
        return len(graph)

    # Create the graph G' = G + v*
    g_prime = add_universal_vertex(graph)

    # We compute C(G, i), the number of i-cliques in G, iteratively.
    # Base case: C(G, 1) is the number of vertices.
    num_cliques_prev = len(graph)
    
    print("Demonstrating the reduction from #k-Clique to PCount.")
    print(f"Base case: The number of 1-cliques is C(G, 1) = {num_cliques_prev}")

    s_values = {}
    for i in range(2, k + 1):
        # Call the oracle for PCount(G', i) to get S_i.
        # For k>=2, PCount(G', i) = C(G, i) + C(G, i-1).
        s_i = p_count_oracle(g_prime, i)
        s_values[i] = s_i
        
        # From C(G, i) + C(G, i-1) = S_i, we get C(G, i) = S_i - C(G, i-1).
        num_cliques_curr = s_i - num_cliques_prev
        
        print(f"For i={i}:")
        print(f"  Oracle for PCount(G', {i}) returns S_{i} = {s_i}")
        print(f"  Number of {i-1}-cliques C(G, {i-1}) = {num_cliques_prev}")
        print(f"  Calculated number of {i}-cliques: C(G, {i}) = S_{i} - C(G, {i-1}) = {s_i} - {num_cliques_prev} = {num_cliques_curr}")
        
        num_cliques_prev = num_cliques_curr

    print("\nThe final result can be expressed via a telescoping sum.")
    print(f"The number of {k}-cliques C(G, {k}) is calculated as {num_cliques_prev}.")
    
    # The full equation is C(G, k) = S_k - S_{k-1} + S_{k-2} - ... + (-1)^(k-1) * C(G, 1)
    print("\nFinal Equation with values:")
    final_eq_values = f"C(G, {k}) = "
    result = 0
    term_list = []

    for i in range(k, 1, -1):
        sign = (-1)**(k-i)
        term = sign * s_values[i]
        result += term
        term_list.append(f"{term: d}")

    sign = (-1)**(k-1)
    term = sign * len(graph)
    result += term
    term_list.append(f"{term: d}")
    
    final_eq_values += " + ".join(term_list).replace("+ -", "- ")
    final_eq_values += f" = {result}"
    print(final_eq_values)

    return num_cliques_prev

# Example usage for demonstration purposes.
if __name__ == '__main__':
    # Example graph: A 4-cycle with a chord, forming one 3-clique.
    # G = {0:[1,2], 1:[0,2], 2:[0,1,3], 3:[2]}
    # Clique counts in G: C(G,1)=4, C(G,2)=4, C(G,3)=1, C(G,4)=0
    example_graph = {
        0: {1, 2},
        1: {0, 2},
        2: {0, 1, 3},
        3: {2}
    }
    k_target = 3

    def dummy_p_count_oracle(graph, k):
        # This oracle simulates the expected behavior for G' = G + v*.
        # PCount(G', k) = C(G, k) + C(G, k-1)
        C_G = {1: 4, 2: 4, 3: 1, 4: 0}
        
        if k > 1:
            return C_G.get(k, 0) + C_G.get(k - 1, 0)
        return C_G.get(1, 0)

    print(f"Attempting to find the number of {k_target}-cliques in the example graph.")
    final_count = count_k_cliques_using_oracle(example_graph, k_target, dummy_p_count_oracle)
    
    print("\nThis demonstrates that if PCount were FPT, #k-Clique would also be FPT.")
    print("Since #k-Clique is #W[1]-hard, PCount must also be #W[1]-hard.")
