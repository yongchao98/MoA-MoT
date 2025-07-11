def solve_super_knight_planarity():
    """
    Calculates the supremum of nm for which the (3,2)-super-knight graph on
    an n x m board (n,m >= 4) is potentially planar.

    The planarity is checked using the necessary condition for bipartite graphs: E <= 2V - 4.
    """

    def calculate_edges(n, m):
        """Calculates the number of edges (E) for the super-knight graph on an n x m board."""
        total_degree = 0
        # All 8 possible moves for a (3,2)-super-knight
        moves = [
            (3, 2), (3, -2), (-3, 2), (-3, -2),
            (2, 3), (2, -3), (-2, 3), (-2, -3)
        ]

        for r in range(n):
            for c in range(m):
                degree = 0
                for dr, dc in moves:
                    nr, nc = r + dr, c + dc
                    # Check if the destination square is on the board
                    if 0 <= nr < n and 0 <= nc < m:
                        degree += 1
                total_degree += degree

        # Each edge is counted twice (once from each vertex), so divide by 2
        return total_degree // 2

    max_nm = 0
    best_config = {}

    # We iterate n from 4 upwards. For each n, we find the largest m that works.
    # The search range is chosen to be large enough to find the maximum.
    # The average degree increases as n and m become more similar (more "squarish"),
    # so we expect the maximum nm to be found for smaller n and larger m.
    for n in range(4, 15):
        # We assume m >= n since the board can be rotated.
        for m in range(n, 50):
            V = n * m
            E = calculate_edges(n, m)
            limit = 2 * V - 4

            if E <= limit:
                # This configuration is a candidate for being planar.
                # If its size is the largest so far, we save it.
                if V > max_nm:
                    max_nm = V
                    best_config = {'n': n, 'm': m, 'V': V, 'E': E, 'limit': limit}
            else:
                # For a fixed n, the ratio of edges to vertices tends to grow with m.
                # Once the condition fails, it will continue to fail for larger m.
                # We can stop searching for this n and move to the next.
                break

    print("The analysis is based on the necessary condition for a bipartite graph to be planar: E <= 2V - 4.")
    print(f"\nThe supremum of the value of nm is {best_config['V']}.")
    print(f"The largest board size found to satisfy the condition is {best_config['n']}x{best_config['m']}.")
    print("\nFor this board configuration:")
    print(f"  Total squares (V) = {best_config['V']}")
    print(f"  Total super-knight moves (E) = {best_config['E']}")
    print(f"  Planarity condition check: E <= 2*V - 4")
    # We output each number in the final equation as requested
    print(f"  {best_config['E']} <= 2*{best_config['V']} - 4")
    print(f"  {best_config['E']} <= {best_config['limit']}")

solve_super_knight_planarity()
