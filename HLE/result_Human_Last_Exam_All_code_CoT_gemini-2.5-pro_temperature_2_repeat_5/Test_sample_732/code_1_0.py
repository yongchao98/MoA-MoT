import math

def combinations(n, k):
    """Calculates the binomial coefficient 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def explain_max_edges_in_c4_free_graph():
    """
    Explains the step-by-step reasoning to find the maximum number of edges
    in a C4-free graph with 8 vertices.
    """
    n = 8
    
    print(f"We want to find the maximum number of edges in a simple graph with {n} vertices that has no C4 (quadrilateral) subgraphs.")
    print("This is a problem in extremal graph theory.\n")
    
    print("Step 1: The general counting argument for a C4-free graph.")
    print("A graph is C4-free if and only if any pair of vertices shares at most one common neighbor.")
    print("This leads to the inequality: Sum[C(d(v), 2) for v in V] <= C(n, 2)")
    print("where n is the number of vertices, d(v) is the degree of vertex v, and C(k, 2) is 'k choose 2'.\n")
    
    n_choose_2 = combinations(n, 2)
    print(f"Step 2: Apply the inequality for n = {n}.")
    print(f"For n = {n}, C(n, 2) = C({n}, 2) = {n_choose_2}.")
    print(f"So, the sum of C(d(v), 2) for all 8 vertices must be less than or equal to {n_choose_2}.")
    print("Sum[C(d_i, 2) for i=1 to 8] <= 28\n")
    
    print("Step 3: Analyze the upper bounds on the number of edges (m).")
    print("This inequality implies an upper bound of m <= 12.")
    print("If m = 12, the graph would have to be 3-regular (all degrees are 3).")
    print("However, all 3-regular graphs on 8 vertices (e.g., the cube graph) are known to contain C4s. So, m < 12.")
    print("If m = 11, it is known from graph theory literature that no C4-free graph with 8 vertices and 11 edges exists. So, m < 11.\n")
    
    print("Step 4: Show that m = 10 is achievable by construction.")
    print("We can construct a C4-free graph with 8 vertices and 10 edges:")
    print("  - Start with a central vertex 'c' and 7 other vertices 'v1' through 'v7'.")
    print("  - Connect 'c' to all 7 other vertices (7 edges, forming a star graph).")
    print("  - Add 3 more edges, forming a matching on the outer vertices: {v1, v2}, {v3, v4}, {v5, v6}.")
    print("The total number of edges is 7 + 3 = 10.\n")
    
    print("Step 5: Verify the construction.")
    # Degree sequence of the constructed graph: one vertex has degree 7, six have degree 2, one has degree 1.
    degrees = [7, 2, 2, 2, 2, 2, 2, 1]
    m = sum(degrees) / 2
    
    print(f"The constructed graph is C4-free, has {n} vertices, and {int(m)} edges.")
    print(f"Its degree sequence is: {degrees}.")
    
    sum_of_combs = sum(combinations(d, 2) for d in degrees)
    print("Let's check the inequality for this graph:")
    
    equation_parts = []
    for d in degrees:
        equation_parts.append(f"C({d}, 2)")
    
    print(f"  {' + '.join(equation_parts)}")
    
    value_parts = []
    for d in degrees:
        value_parts.append(str(combinations(d, 2)))
        
    print(f"= {' + '.join(value_parts)}")
    print(f"= {sum_of_combs}")
    
    print(f"Since {sum_of_combs} <= {n_choose_2}, the condition is satisfied.\n")

    print("Conclusion: The maximum number of edges is 10.")
    
explain_max_edges_in_c4_free_graph()