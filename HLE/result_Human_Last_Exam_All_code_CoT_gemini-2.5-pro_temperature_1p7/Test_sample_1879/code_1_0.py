import collections

def solve_grid_paths():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_TARGET = 4
    U_TARGET = 8
    MAX_CONSECUTIVE = 3

    # dp[r][u][dir][k]: num ways to get to (r, u)
    # ending with k consecutive moves in direction dir (0=Up, 1=Right).
    dp = [[[[0] * (MAX_CONSECUTIVE + 1) for _ in range(2)] for _ in range(U_TARGET + 1)] for _ in range(R_TARGET + 1)]

    # Base case: At (0,0), there is one path of length 0.
    # We imagine this path can transition to the first step.
    # To simplify, we'll set the states for (1,0) and (0,1) directly.
    if R_TARGET > 0:
        dp[1][0][1][1] = 1 # Path "R"
    if U_TARGET > 0:
        dp[0][1][0][1] = 1 # Path "U"

    # Iterate through the grid, skipping (0,0) and the base cases
    for r in range(R_TARGET + 1):
        for u in range(U_TARGET + 1):
            if (r == 0 and u == 0) or (r == 1 and u == 0) or (r == 0 and u == 1):
                continue
            
            # Calculate paths ending with an UP move (from (r, u-1))
            if u > 0:
                # To end with 1 UP move, path to (r, u-1) must have ended in RIGHT moves.
                paths_from_R = sum(dp[r][u-1][1])
                dp[r][u][0][1] = paths_from_R
                
                # To end with k > 1 UP moves, path to (r, u-1) must have ended in k-1 UP moves.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    dp[r][u][0][k] = dp[r][u-1][0][k-1]

            # Calculate paths ending with a RIGHT move (from (r-1, u))
            if r > 0:
                # To end with 1 RIGHT move, path to (r-1, u) must have ended in UP moves.
                paths_from_U = sum(dp[r-1][u][0])
                dp[r][u][1][1] = paths_from_U
                
                # To end with k > 1 RIGHT moves, path to (r-1, u) must have ended in k-1 RIGHT moves.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    dp[r][u][1][k] = dp[r-1][u][1][k-1]

    # Get the results for the target cell (4,8)
    final_U_paths = dp[R_TARGET][U_TARGET][0]
    final_R_paths = dp[R_TARGET][U_TARGET][1]
    
    total_U = sum(final_U_paths)
    total_R = sum(final_R_paths)
    
    print("Paths to (4,8) can end in one of two ways:")
    print(f"1. Ending with one or more UP moves: {final_U_paths[1]} (U) + {final_U_paths[2]} (UU) + {final_U_paths[3]} (UUU) = {total_U}")
    print(f"2. Ending with one or more RIGHT moves: {final_R_paths[1]} (R) + {final_R_paths[2]} (RR) + {final_R_paths[3]} (RRR) = {total_R}")
    
    total_ways = total_U + total_R
    print(f"\nTotal unique ways = {total_U} + {total_R} = {total_ways}")

solve_grid_paths()
<<<155>>>