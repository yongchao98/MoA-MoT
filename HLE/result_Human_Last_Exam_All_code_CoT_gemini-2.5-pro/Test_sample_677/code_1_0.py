import sys

def count_edges(n, m):
    """Counts the number of edges in the (3,2) super-knight graph on an n x m board."""
    edges = 0
    # All 8 possible (3,2) super-knight moves
    moves = [(3, 2), (3, -2), (-3, 2), (-3, -2),
             (2, 3), (2, -3), (-2, 3), (-2, -3)]
    
    # Since each edge is counted twice (once from each endpoint), we sum degrees and divide by 2.
    total_degree = 0
    for r in range(n):
        for c in range(m):
            degree = 0
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                # Check if the destination square is on the board
                if 0 <= nr < n and 0 <= nc < m:
                    degree += 1
            total_degree += degree
            
    return total_degree // 2

def find_non_planar_boundary():
    """
    Finds the smallest board areas (nm) with n,m>=4 that violate the
    necessary planarity condition for bipartite graphs (E <= 2V - 4).
    """
    print("Searching for non-planar boards (n, m >= 4)...")
    print("A board is non-planar if E > 2V - 4.\n")
    
    violations = []
    # Search a reasonable range of n and m
    for n in range(4, 15):
        for m in range(n, 20):
            V = n * m
            E = count_edges(n, m)
            limit = 2 * V - 4
            
            if E > limit:
                violations.append({'n': n, 'm': m, 'nm': V, 'E': E, 'limit': limit})

    # Sort violations by area nm to find the smallest ones
    sorted_violations = sorted(violations, key=lambda x: x['nm'])
    
    print("Found the following boards that are provably non-planar:")
    for v in sorted_violations[:5]: # Print the top 5 smallest
        print(f"Board: {v['n']}x{v['m']}, Area nm={v['nm']}. "
              f"Equation: E = {v['E']}, 2V - 4 = {v['limit']}. "
              f"Result: {v['E']} > {v['limit']}")

find_non_planar_boundary()
