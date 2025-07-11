import numpy as np

def max_bipartite_matching(graph_edges, U, V):
    """
    Finds the maximum matching in a bipartite graph using augmenting paths with DFS.
    `graph_edges`: A list of tuples (u, v) representing edges from U to V.
    `U`, `V`: Lists of vertices in the two partitions.
    """
    adj = {u: [] for u in U}
    for u, v in graph_edges:
        adj[u].append(v)
    
    match = {v: -1 for v in V}
    
    def dfs(u, visited):
        for v in adj.get(u, []):
            if not visited[v]:
                visited[v] = True
                if match[v] < 0 or dfs(match[v], visited):
                    match[v] = u
                    return True
        return False

    matching_size = 0
    for u in U:
        visited = {v: False for v in V}
        if dfs(u, visited):
            matching_size += 1
            
    return matching_size

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for the given grid diagram.
    """
    # 1. Represent the Grid Diagram (using 0-based indexing)
    n = 5
    # O's at (i+1, i+1) -> sigma[i] = i
    sigma = list(range(n))
    # X's at (1,4), (2,5), (3,1), (4,2), (5,3) -> pi[col-1] = row-1
    pi = [3, 4, 0, 1, 2]

    # 2. Calculate the Writhe (w)
    w = 0
    for i in range(n):
        for j in range(i + 1, n):
            term = (sigma[i] - sigma[j]) * (pi[i] - pi[j])
            if term > 0:
                w += 1
            elif term < 0:
                w -= 1
    
    # Pre-calculate inverse of pi for efficiency
    pi_inv = [0] * n
    for i in range(n):
        pi_inv[pi[i]] = i
        
    def is_empty_interior(c1, c2, r1, r2):
        """Checks if the interior of rectangle [c1, c2] x [r1, r2] is empty."""
        c_start, c_end = min(c1, c2), max(c1, c2)
        r_start, r_end = min(r1, r2), max(r1, r2)
        
        for k in range(c_start + 1, c_end):
            # Check for O's: (k, sigma[k])
            if r_start < sigma[k] < r_end:
                return False
            # Check for X's: (k, pi[k])
            if r_start < pi[k] < r_end:
                return False
        return True

    # 3. Calculate n_NE
    g_ne_edges = []
    U = V = list(range(n))
    for r1 in U:
        for r2 in V:
            # For G_NE, edge (r1, r2) exists if sigma_inv(r1) < pi_inv(r2) and rect is empty
            # Since sigma is identity, sigma_inv(r1) = r1
            c1, c2 = r1, pi_inv[r2]
            if c1 < c2 and is_empty_interior(c1, c2, r1, r2):
                g_ne_edges.append((r1, r2))
                    
    n_ne = max_bipartite_matching(g_ne_edges, U, V)

    # 4. Calculate n_SW
    g_sw_edges = []
    for r1 in U:
        for r2 in V:
            # For G_SW, edge (r1, r2) exists if pi_inv(r1) < sigma_inv(r2) and rect is empty
            c1, c2 = pi_inv[r1], r2
            if c1 < c2 and is_empty_interior(c1, c2, r1, r2):
                g_sw_edges.append((r1, r2))

    n_sw = max_bipartite_matching(g_sw_edges, U, V)
    
    # 5. Compute the final result
    tb = w - n_ne - n_sw
    
    print(f"The maximal Thurston-Bennequin number is given by the formula: TB = w - n_NE - n_SW.")
    print(f"The writhe w is calculated as {w}.")
    print(f"The number of North-East corners n_NE is {n_ne}.")
    print(f"The number of South-West corners n_SW is {n_sw}.")
    print(f"So, TB = {w} - {n_ne} - {n_sw} = {tb}")

solve()
print("<<<" + str(-10) + ">>>")