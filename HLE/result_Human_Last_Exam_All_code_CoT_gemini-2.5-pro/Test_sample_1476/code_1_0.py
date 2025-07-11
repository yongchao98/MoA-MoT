import numpy as np

def calculate_total_variation():
    """
    This function demonstrates the concept based on the problem's deductions.
    The reasoning implies that the edge signal x1 must be the zero vector.
    This happens if and only if the vertex signal x0 is constant across
    all connected vertices, which means the graph has a total variation of 0.
    """
    # Let's define a graph G as a set of edges.
    # Consider a square graph with 4 vertices (0, 1, 2, 3) and 4 edges.
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    num_vertices = 4

    # The conclusion x1 = 0 implies |x0_u - x0_v| = 0 for all edges {u,v}.
    # This means x0 must be constant on connected components of the graph.
    # Let's define such a vertex signal x0.
    x0 = [12.5, 12.5, 12.5, 12.5]
    print(f"Let G be a graph with vertices V and edges E.")
    print(f"Let the vertex signal be x0 = {x0}\n")
    print(f"The total variation of x0 on G is defined as:")
    print(f"TV(x0) = sum_{{u,v} in E} |x0_u - x0_v|\n")

    # Now, let's calculate the total variation by summing the absolute differences
    # for each edge, and print the equation.
    
    equation_terms = []
    value_terms = []
    total_variation = 0.0

    for u, v in edges:
        diff = np.abs(x0[u] - x0[v])
        equation_terms.append(f"|{x0[u]} - {x0[v]}|")
        value_terms.append(str(diff))
        total_variation += diff

    # Print the equation with all the numbers
    equation_str = " + ".join(equation_terms)
    values_str = " + ".join(value_terms)
    
    print("For our graph and signal, the calculation is:")
    print(f"TV(x0) = {equation_str}")
    print(f"       = {values_str}")
    print(f"       = {total_variation}")

    print("\nAs we can see, the total variation is 0.")
    print("This confirms the inference from the problem statement.")

# Run the calculation and print the results
calculate_total_variation()