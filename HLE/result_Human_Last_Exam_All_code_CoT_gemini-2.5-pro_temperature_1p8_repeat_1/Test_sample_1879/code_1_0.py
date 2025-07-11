def solve_grid_paths():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with constraints on consecutive moves using dynamic programming.
    """
    R_target = 4
    U_target = 8
    K_max = 3  # Maximum allowed consecutive moves in one direction

    # DP tables
    # R_ways[r][u][k]: number of paths to (r,u) ending in k consecutive Right moves.
    # U_ways[r][u][k]: number of paths to (r,u) ending in k consecutive Up moves.
    # k is 1-indexed (k=1, 2, 3), so the array size is K_max + 1.
    R_ways = [[[0] * (K_max + 1) for _ in range(U_target + 1)] for _ in range(R_target + 1)]
    U_ways = [[[0] * (K_max + 1) for _ in range(U_target + 1)] for _ in range(R_target + 1)]

    # Loop to fill the DP tables for each point (r, u) on the grid.
    for r in range(R_target + 1):
        for u in range(U_target + 1):
            if r == 0 and u == 0:
                continue

            # Calculate ways for paths ending with 'Right' moves to reach (r, u)
            if r > 0:
                # Paths ending with a single 'R' must come from a path ending in 'U'.
                # Sum all ways to reach (r-1, u) that ended in any 'U' sequence.
                sum_from_up = sum(U_ways[r - 1][u])
                
                # Base Case: The very first step 'R' from the origin (0,0) to (1,0).
                # This is a special case as (0,0) has no preceding move.
                # We model this as 1 way to start with 'R'.
                if r == 1 and u == 0:
                    sum_from_up = 1
                    
                R_ways[r][u][1] = sum_from_up
                
                # Paths ending with k > 1 'R's must extend a path that ended in k-1 'R's.
                for k in range(2, K_max + 1):
                    R_ways[r][u][k] = R_ways[r - 1][u][k - 1]

            # Calculate ways for paths ending with 'Up' moves to reach (r, u)
            if u > 0:
                # Paths ending with a single 'U' must come from a path ending in 'R'.
                # Sum all ways to reach (r, u-1) that ended in any 'R' sequence.
                sum_from_right = sum(R_ways[r][u - 1])

                # Base Case: The very first step 'U' from the origin (0,0) to (0,1).
                if r == 0 and u == 1:
                    sum_from_right = 1

                U_ways[r][u][1] = sum_from_right

                # Paths ending with k > 1 'U's must extend a path that ended in k-1 'U's.
                for k in range(2, K_max + 1):
                    U_ways[r][u][k] = U_ways[r][u - 1][k - 1]

    # The final answer is the sum of all valid paths reaching the target B(4,8).
    final_R_ways = R_ways[R_target][U_target]
    final_U_ways = U_ways[R_target][U_target]
    
    total_ways = sum(final_R_ways) + sum(final_U_ways)

    print("The total number of unique ways is the sum of paths ending in a Right move and paths ending in an Up move.")
    print("\nCalculation for paths to B(4,8):")
    print("-" * 35)

    r_comp_str = " + ".join(map(str, final_R_ways[1:]))
    print(f"Paths ending in Right moves: {r_comp_str} = {sum(final_R_ways)}")
    
    u_comp_str = " + ".join(map(str, final_U_ways[1:]))
    print(f"Paths ending in Up moves:    {u_comp_str} = {sum(final_U_ways)}")
    
    print("-" * 35)
    print(f"Total Unique Ways = {sum(final_R_ways)} + {sum(final_U_ways)} = {total_ways}")

solve_grid_paths()
<<<79>>>