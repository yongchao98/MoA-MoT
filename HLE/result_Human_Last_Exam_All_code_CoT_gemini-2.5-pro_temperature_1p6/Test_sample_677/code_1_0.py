import sys

def calculate_edges(n, m):
    """
    Calculates the number of edges in the super-knight graph on an n x m board.
    """
    moves = [
        (2, 3), (2, -3), (-2, 3), (-2, -3),
        (3, 2), (3, -2), (-3, 2), (-3, -2)
    ]
    degree_sum = 0
    for r in range(n):
        for c in range(m):
            degree = 0
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < m:
                    degree += 1
            degree_sum += degree
    # Each edge is counted twice (once from each endpoint vertex)
    return degree_sum // 2

def find_supremum_size():
    """
    Finds the maximum board size nm that satisfies the planar graph condition E <= 2V - 4.
    """
    max_size = 0
    n_at_max = 0
    m_at_max = 0
    E_at_max = 0
    V_at_max = 0
    
    # We search up to a reasonable limit for n and m.
    # If a board (n, m) becomes non-planar by this criterion, (n, m+1) will also be non-planar,
    # as the number of edges grows faster than the number of vertices.
    # This allows us to break the inner loop early.
    limit = 40  # A sufficiently large limit for the search
    for n in range(4, limit):
        for m in range(n, limit):
            V = n * m
            E = calculate_edges(n, m)
            
            # Planarity condition for bipartite graphs
            if E <= 2 * V - 4:
                if V > max_size:
                    max_size = V
                    n_at_max = n
                    m_at_max = m
                    E_at_max = E
                    V_at_max = V
            else:
                # Break inner loop as subsequent m for this n will also fail
                break
                
    print(f"The supremum of the value of the largest size found is {max_size}.")
    print(f"This occurs for a {n_at_max}x{m_at_max} board (and its {m_at_max}x{n_at_max} transpose).")
    print("For this board, we check the planarity condition E <= 2V - 4.")
    print("The final equation is:")
    print(f"{E_at_max} <= 2 * {V_at_max} - 4")


find_supremum_size()
