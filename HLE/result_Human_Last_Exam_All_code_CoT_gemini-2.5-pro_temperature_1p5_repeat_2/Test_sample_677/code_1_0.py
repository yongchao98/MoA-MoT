def solve_superknight_planarity():
    """
    Determines the supremum of the area nm for an n x m chessboard (n, m >= 4)
    where the (3,2) super-knight graph is planar.
    
    The method is to find the largest area nm that satisfies the necessary planarity
    condition for graphs without C3 cycles: e <= 2v - 4.
    """
    
    # All 8 possible super-knight moves
    moves = [
        (3, 2), (3, -2), (-3, 2), (-3, -2),
        (2, 3), (2, -3), (-2, 3), (-2, -3)
    ]

    def count_edges(n, m):
        """Counts the total number of edges in the super-knight graph for an n x m board."""
        total_degree = 0
        for r in range(n):
            for c in range(m):
                degree = 0
                for dr, dc in moves:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < n and 0 <= nc < m:
                        degree += 1
                total_degree += degree
        # Each edge is counted twice (once from each vertex), so we divide by 2.
        return total_degree // 2

    max_planar_area = 0
    best_n, best_m = 0, 0
    
    # Search up to a reasonable limit for n and m
    search_limit_n = 15
    search_limit_m = 30
    
    print("Searching for the largest area 'nm' satisfying the planarity condition e <= 2v - 4...")
    print("-" * 70)

    for n in range(4, search_limit_n):
        for m in range(n, search_limit_m):
            v = n * m
            e = count_edges(n, m)
            
            # The necessary condition for a planar graph with no triangles
            planar_condition_holds = (e <= 2 * v - 4)
            
            if planar_condition_holds:
                if v > max_planar_area:
                    max_planar_area = v
                    best_n, best_m = n, m
            else:
                # Since e grows faster than 2v, larger boards for this 'n' will also fail.
                # We can break the inner loop and proceed to the next 'n'.
                break
                
    print(f"\nSearch complete.")
    print(f"The maximum area found that satisfies the condition is {max_planar_area}.")
    print(f"This occurs for a board of size {best_n} x {best_m} (or {best_m} x {best_n}).")
    
    # Finally, show the calculation for the optimal board size found
    final_n, final_m = best_n, best_m
    final_v = final_n * final_m
    final_e = count_edges(final_n, final_m)
    limit = 2 * final_v - 4
    
    print("\nVerification for the board with the maximum planar area:")
    print(f"Board size: n = {final_n}, m = {final_m}")
    print(f"Number of vertices (v): {final_v}")
    print(f"Number of edges (e): {final_e}")
    print("\nChecking the planarity condition e <= 2v - 4:")
    print(f"{final_e} <= 2 * {final_v} - 4")
    print(f"{final_e} <= {limit}")
    print(f"The condition holds.")

solve_superknight_planarity()