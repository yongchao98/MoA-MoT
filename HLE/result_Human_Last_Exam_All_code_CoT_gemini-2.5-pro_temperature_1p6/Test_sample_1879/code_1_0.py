def solve_grid_paths():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with specific movement constraints using dynamic programming.
    """
    target_r = 4
    target_u = 8
    max_consecutive = 3

    # dp[i][j][d][c]:
    # i: number of right steps
    # j: number of up steps
    # d: direction of the last move (0 for Right, 1 for Up)
    # c: consecutive steps in that direction (1, 2, or 3)
    # Initialize a 4D list with zeros
    dp = [[[([0] * (max_consecutive + 1)) for _ in range(2)] for _ in range(target_u + 1)] for _ in range(target_r + 1)]

    # Iterate through each cell of the grid
    for i in range(target_r + 1):
        for j in range(target_u + 1):
            if i == 0 and j == 0:
                continue

            # --- Calculate ways ending in RIGHT moves ---
            # A path to (i,j) ending in R must come from (i-1, j)
            if i > 0:
                # Case 1: The path ends in a single R (e.g., ...UR)
                # This requires the path to (i-1,j) to have ended in U.
                # We sum up all ways to (i-1,j) that ended in 1, 2, or 3 U's.
                sum_from_u = 0
                if i == 1 and j == 0:
                    # This is the very first step being 'R' from the start.
                    sum_from_u = 1
                else:
                    for k in range(1, max_consecutive + 1):
                        sum_from_u += dp[i-1][j][1][k]
                dp[i][j][0][1] = sum_from_u

                # Case 2 & 3: The path ends in RR or RRR.
                # A path ending in c consecutive R's must come from a path
                # that ended in c-1 consecutive R's at the previous step.
                for c in range(2, max_consecutive + 1):
                    dp[i][j][0][c] = dp[i-1][j][0][c-1]
            
            # --- Calculate ways ending in UP moves ---
            # A path to (i,j) ending in U must come from (i, j-1)
            if j > 0:
                # Case 1: The path ends in a single U (e.g., ...RU)
                # This requires the path to (i, j-1) to have ended in R.
                sum_from_r = 0
                if i == 0 and j == 1:
                    # This is the very first step being 'U' from the start.
                    sum_from_r = 1
                else:
                    for k in range(1, max_consecutive + 1):
                        sum_from_r += dp[i][j-1][0][k]
                dp[i][j][1][1] = sum_from_r
            
                # Case 2 & 3: The path ends in UU or UUU.
                for c in range(2, max_consecutive + 1):
                    dp[i][j][1][c] = dp[i][j-1][1][c-1]

    # Extract the final results for B(4,8)
    # Number of ways ending in 1, 2, or 3 Right moves
    r1 = dp[target_r][target_u][0][1]
    r2 = dp[target_r][target_u][0][2]
    r3 = dp[target_r][target_u][0][3]
    
    # Number of ways ending in 1, 2, or 3 Up moves
    u1 = dp[target_r][target_u][1][1]
    u2 = dp[target_r][target_u][1][2]
    u3 = dp[target_r][target_u][1][3]

    total_ways = r1 + r2 + r3 + u1 + u2 + u3

    print("To find the total number of unique ways, we sum the ways to arrive at (4,8) ending with each possible valid sequence of moves:")
    print(f"Ways ending in ...U(R): {r1}")
    print(f"Ways ending in ..UR(R): {r2}")
    print(f"Ways ending in .URR(R): {r3}")
    print(f"Ways ending in ...R(U): {u1}")
    print(f"Ways ending in ..RU(U): {u2}")
    print(f"Ways ending in .RUU(U): {u3}")
    print("-" * 20)
    print(f"Final calculation: {r1} + {r2} + {r3} + {u1} + {u2} + {u3} = {total_ways}")

solve_grid_paths()