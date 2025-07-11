def solve_grid_path_problem():
    """
    Calculates the number of unique paths from (0,0) to (R,U) on a grid with constraints.
    Movement is restricted to 1 unit right or 1 unit up per step.
    A sequence of four or more consecutive steps in the same direction is forbidden.
    """
    # Grid dimensions and constraints
    R = 4  # Target x-coordinate (number of Right moves)
    U = 8  # Target y-coordinate (number of Up moves)
    K = 3  # Maximum allowed consecutive steps in the same direction

    # DP table: dp[r][u][dir][k]
    # r: number of Right moves (0 to R)
    # u: number of Up moves (0 to U)
    # dir: 0 for Right, 1 for Up
    # k: number of consecutive steps in 'dir' (1-based index: 1, 2, 3)
    dp = [[[[0] * (K + 1) for _ in range(2)] for _ in range(U + 1)] for _ in range(R + 1)]

    # Helper table: ways[r][u][dir]
    # Stores the total number of ways to reach (r,u) ending with any number
    # of consecutive moves in direction 'dir'.
    ways = [[[0] * 2 for _ in range(U + 1)] for _ in range(R + 1)]

    # Seed the calculation at the origin (0,0).
    # This is a conceptual trick to make the DP transitions work from the start.
    # It represents the choice of making the first step Right or Up from a virtual start point.
    ways[0][0][0] = 1
    ways[0][0][1] = 1

    # Iterate over the grid to fill the DP tables
    for r in range(R + 1):
        for u in range(U + 1):
            if r == 0 and u == 0:
                continue

            # Calculate ways to reach (r, u) ending with a Right move
            if r > 0:
                # To end with 1 Right move, the previous path must have ended with an Up move.
                dp[r][u][0][1] = ways[r-1][u][1]
                # To end with k > 1 Right moves, the previous path must have ended with k-1 Right moves.
                for k in range(2, K + 1):
                    dp[r][u][0][k] = dp[r-1][u][0][k-1]
                # Update the total ways ending in Right for the current cell
                ways[r][u][0] = sum(dp[r][u][0])

            # Calculate ways to reach (r, u) ending with an Up move
            if u > 0:
                # To end with 1 Up move, the previous path must have ended with a Right move.
                dp[r][u][1][1] = ways[r][u-1][0]
                # To end with k > 1 Up moves, the previous path must have ended with k-1 Up moves.
                for k in range(2, K + 1):
                    dp[r][u][1][k] = dp[r][u-1][1][k-1]
                # Update the total ways ending in Up for the current cell
                ways[r][u][1] = sum(dp[r][u][1])

    # The total number of ways is the sum of paths ending in Right or Up at the destination (R,U)
    total_ways = ways[R][U][0] + ways[R][U][1]

    # Breakdown of the final sum
    w_r_1 = dp[R][U][0][1]
    w_r_2 = dp[R][U][0][2]
    w_r_3 = dp[R][U][0][3]
    w_u_1 = dp[R][U][1][1]
    w_u_2 = dp[R][U][1][2]
    w_u_3 = dp[R][U][1][3]

    print(f"The number of unique ways from (0,0) to ({R},{U}) is found by summing the ways to arrive with different final steps:")
    print(f"Ways ending in 1 Right step: {w_r_1}")
    print(f"Ways ending in 2 Right steps: {w_r_2}")
    print(f"Ways ending in 3 Right steps: {w_r_3}")
    print(f"Ways ending in 1 Up step: {w_u_1}")
    print(f"Ways ending in 2 Up steps: {w_u_2}")
    print(f"Ways ending in 3 Up steps: {w_u_3}")
    print(f"\nFinal Equation:")
    print(f"Total ways = {w_r_1} + {w_r_2} + {w_r_3} + {w_u_1} + {w_u_2} + {w_u_3} = {total_ways}")
    print(f"\n<<< {total_ways} >>>")

solve_grid_path_problem()