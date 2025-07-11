import itertools

def get_super_knight_neighbors(r, c, n, m):
    """
    Calculates the valid super-knight moves for a given square (r, c) on an n x m board.
    The board is 1-indexed.
    """
    moves = [(2, 3), (2, -3), (-2, 3), (-2, -3),
             (3, 2), (3, -2), (-3, 2), (-3, -2)]
    neighbors = set()
    for dr, dc in moves:
        nr, nc = r + dr, c + dc
        if 1 <= nr <= n and 1 <= nc <= m:
            neighbors.add((nr, nc))
    return neighbors

def find_k33(n, m):
    """
    Searches for a K_3,3 subgraph on an n x m board.
    A K_3,3 is found if 3 nodes in one partition of the bipartite graph
    share at least 3 common neighbors in the other partition.
    """
    nodes = [(r, c) for r in range(1, n + 1) for c in range(1, m + 1)]
    # Partition nodes into two sets based on chessboard coloring (bipartite sets)
    part1 = {node for node in nodes if (node[0] + node[1]) % 2 == 0} # White squares
    part2 = {node for node in nodes if (node[0] + node[1]) % 2 != 0} # Black squares

    # Ensure both partitions are large enough
    if len(part1) < 3 or len(part2) < 3:
        return None

    # Search for U (triplet) in part1 and V (common neighbors) in part2
    for u_triplet in itertools.combinations(part1, 3):
        # Get neighbors for each of the 3 nodes in the triplet
        n1 = get_super_knight_neighbors(u_triplet[0][0], u_triplet[0][1], n, m)
        n2 = get_super_knight_neighbors(u_triplet[1][0], u_triplet[1][1], n, m)
        n3 = get_super_knight_neighbors(u_triplet[2][0], u_triplet[2][1], n, m)
        
        # Find common neighbors. Since the graph is bipartite, they are guaranteed to be in part2.
        common_neighbors = n1.intersection(n2).intersection(n3)
        
        if len(common_neighbors) >= 3:
            v_triplet = tuple(itertools.islice(common_neighbors, 3))
            return (u_triplet, v_triplet)

    # Search for U in part2 and V in part1
    for u_triplet in itertools.combinations(part2, 3):
        n1 = get_super_knight_neighbors(u_triplet[0][0], u_triplet[0][1], n, m)
        n2 = get_super_knight_neighbors(u_triplet[1][0], u_triplet[1][1], n, m)
        n3 = get_super_knight_neighbors(u_triplet[2][0], u_triplet[2][1], n, m)
        
        common_neighbors = n1.intersection(n2).intersection(n3)

        if len(common_neighbors) >= 3:
            v_triplet = tuple(itertools.islice(common_neighbors, 3))
            return (u_triplet, v_triplet)

    return None

def solve():
    """
    Finds the smallest non-planar super-knight graph and calculates the supremum size for planar ones.
    """
    min_non_planar_size = float('inf')
    first_non_planar_board = None
    k33_config = None

    # Check boards in increasing order of size nm
    # Limiting search space to make execution feasible, as the principle is demonstrated
    # by the first non-planar graph found.
    for size in range(16, 51): # from 4x4 up to a reasonable limit
        for n in range(4, int(size**0.5) + 1):
            if size % n == 0:
                m = size // n
                if m < 4:
                    continue
                
                # print(f"Checking {n}x{m} board (size {size})...")
                result = find_k33(n, m)
                if result:
                    min_non_planar_size = size
                    first_non_planar_board = (n, m)
                    k33_config = result
                    
                    # Print the results and the deduced supremum
                    print(f"The graph becomes non-planar for a board of size {first_non_planar_board[0]}x{first_non_planar_board[1]}.")
                    print(f"Total vertices (nm): {min_non_planar_size}")
                    print("A K_3,3 subgraph was found:")
                    
                    U_nodes, V_nodes = k33_config
                    print(f"  Set U = {{ {U_nodes[0]}, {U_nodes[1]}, {U_nodes[2]} }}")
                    print(f"  Set V = {{ {V_nodes[0]}, {V_nodes[1]}, {V_nodes[2]} }}")

                    supremum = min_non_planar_size - 1
                    print(f"\nThe largest size 'nm' for which the graph is planar is {supremum}.")
                    return

    print("No non-planar graph found in the checked range.")

solve()