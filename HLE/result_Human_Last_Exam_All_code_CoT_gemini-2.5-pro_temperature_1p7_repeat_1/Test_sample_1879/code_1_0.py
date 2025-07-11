def solve_grid_walk():
    """
    Calculates the number of unique paths on a 2D grid from (0,0) to (4,8)
    with constraints on movement.
    """
    # Grid dimensions and constraints
    R_MAX = 4
    U_MAX = 8
    CONSECUTIVE_MAX = 3

    # DP table definition: dp[r][u][dir][k]
    # r: right steps (x-coordinate), from 0 to R_MAX
    # u: up steps (y-coordinate), from 0 to U_MAX
    # dir: last direction (0 for Right, 1 for Up)
    # k: consecutive moves in last direction (1 to CONSECUTIVE_MAX)
    # We use k from 1 to 3, so the size of the last dimension is 4 (index 0 is unused).
    dp = [[[[0] * (CONSECUTIVE_MAX + 1) for _ in range(2)] for _ in range(U_MAX + 1)] for _ in range(R_MAX + 1)]

    # Iterate through each cell of the grid to compute the number of ways to reach it.
    for r in range(R_MAX + 1):
        for u in range(U_MAX + 1):
            if r == 0 and u == 0:
                continue

            # --- Calculate ways to arrive at (r, u) with a RIGHT move ---
            if r > 0:
                # Case 1: Ending with k=1 right move (e.g., ...UR)
                # This requires the path to (r-1, u) to end with one or more UP moves.
                # A special case is the very first move from (0,0) to (1,0).
                if r == 1 and u == 0:
                    ways_from_up = 1  # The path is just "R"
                else:
                    ways_from_up = sum(dp[r - 1][u][1])
                dp[r][u][0][1] = ways_from_up

                # Case 2: Ending with k > 1 right moves (e.g., ...RR, ...RRR)
                # To end with k 'R's at (r, u), path to (r-1, u) must have ended with k-1 'R's.
                for k in range(2, CONSECUTIVE_MAX + 1):
                    dp[r][u][0][k] = dp[r - 1][u][0][k - 1]

            # --- Calculate ways to arrive at (r, u) with an UP move ---
            if u > 0:
                # Case 1: Ending with k=1 up move (e.g., ...RU)
                # This requires the path to (r, u-1) to end with one or more RIGHT moves.
                # A special case is the very first move from (0,0) to (0,1).
                if r == 0 and u == 1:
                    ways_from_right = 1  # The path is just "U"
                else:
                    ways_from_right = sum(dp[r][u - 1][0])
                dp[r][u][1][1] = ways_from_right
                
                # Case 2: Ending with k > 1 up moves (e.g., ...UU, ...UUU)
                # To end with k 'U's at (r, u), path to (r, u-1) must have ended with k-1 'U's.
                for k in range(2, CONSECUTIVE_MAX + 1):
                    dp[r][u][1][k] = dp[r][u - 1][1][k - 1]

    # The total number of ways is the sum of all possible states at the destination B(4,8).
    final_ways_ending_R = dp[R_MAX][U_MAX][0]
    final_ways_ending_U = dp[R_MAX][U_MAX][1]
    
    total_ways = sum(final_ways_ending_R) + sum(final_ways_ending_U)

    # Output the detailed breakdown of the final answer as an equation.
    print("The total number of unique ways is the sum of ways to arrive at B(4,8) ending with any valid sequence of moves.")
    print("Paths can end with:")
    print(f"- 1, 2, or 3 consecutive RIGHT moves: {final_ways_ending_R[1]}, {final_ways_ending_R[2]}, and {final_ways_ending_R[3]} ways respectively.")
    print(f"- 1, 2, or 3 consecutive UP moves: {final_ways_ending_U[1]}, {final_ways_ending_U[2]}, and {final_ways_ending_U[3]} ways respectively.")
    print("\nThe final equation is:")
    print(f"{final_ways_ending_R[1]} + {final_ways_ending_R[2]} + {final_ways_ending_R[3]} + {final_ways_ending_U[1]} + {final_ways_ending_U[2]} + {final_ways_ending_U[3]} = {total_ways}")
    
    return total_ways

if __name__ == '__main__':
    solve_grid_walk()