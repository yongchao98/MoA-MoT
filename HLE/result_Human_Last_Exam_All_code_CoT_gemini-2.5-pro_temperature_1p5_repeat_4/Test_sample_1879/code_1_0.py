import sys

def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (R_max, U_max)
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_max = 4
    U_max = 8

    # dp[r][u][0..2] will store paths ending in 1, 2, or 3 'R' moves (h1, h2, h3).
    # dp[r][u][3..5] will store paths ending in 1, 2, or 3 'U' moves (v1, v2, v3).
    # We use a list of 6 integers for the last dimension.
    dp = [[[0] * 6 for _ in range(U_max + 1)] for _ in range(R_max + 1)]

    # Base case: Path 'R' to (1,0). This path ends in one 'R'.
    if R_max >= 1:
        dp[1][0][0] = 1

    # Base case: Path 'U' to (0,1). This path ends in one 'U'.
    if U_max >= 1:
        dp[0][1][3] = 1

    # Fill the DP table by iterating through all grid cells.
    for r in range(R_max + 1):
        for u in range(U_max + 1):
            # Skip origin and base cases which are already set.
            if (r == 0 and u == 0) or (r == 1 and u == 0) or (r == 0 and u == 1):
                continue
            
            # Calculate paths ending in 'R' moves for cell (r,u).
            # These must come from cell (r-1, u).
            if r > 0:
                # To end in 1 'R', the previous path at (r-1,u) must have ended in 'U's.
                total_from_U = dp[r-1][u][3] + dp[r-1][u][4] + dp[r-1][u][5]
                dp[r][u][0] = total_from_U
                # To end in 2 'R's, the previous path at (r-1,u) must have ended in 1 'R'.
                dp[r][u][1] = dp[r-1][u][0]
                # To end in 3 'R's, the previous path at (r-1,u) must have ended in 2 'R's.
                dp[r][u][2] = dp[r-1][u][1]

            # Calculate paths ending in 'U' moves for cell (r,u).
            # These must come from cell (r, u-1).
            if u > 0:
                # To end in 1 'U', the previous path at (r,u-1) must have ended in 'R's.
                total_from_R = dp[r][u-1][0] + dp[r][u-1][1] + dp[r][u-1][2]
                dp[r][u][3] = total_from_R
                # To end in 2 'U's, the previous path at (r,u-1) must have ended in 1 'U'.
                dp[r][u][4] = dp[r][u-1][3]
                # To end in 3 'U's, the previous path at (r,u-1) must have ended in 2 'U's.
                dp[r][u][5] = dp[r][u-1][4]

    # The final answer is the sum of all possible ways to arrive at (R_max, U_max).
    final_counts = dp[R_max][U_max]
    total_ways = sum(final_counts)
    
    # Unpack for formatted printing
    h1, h2, h3, v1, v2, v3 = final_counts
    
    print(f"Number of ways to reach ({R_max},{U_max}) ending in:")
    print(f"  - 1 Right move: {h1}")
    print(f"  - 2 Right moves: {h2}")
    print(f"  - 3 Right moves: {h3}")
    print(f"  - 1 Up move: {v1}")
    print(f"  - 2 Up moves: {v2}")
    print(f"  - 3 Up moves: {v3}")
    print("\nTotal number of unique ways:")
    print(f"{h1} + {h2} + {h3} + {v1} + {v2} + {v3} = {total_ways}")
    
    # Final answer for parsing
    # To prevent any extra text from interfering with the final answer extraction,
    # we print it last, with a clear marker.
    # We flush stdout to ensure order of printing is preserved
    sys.stdout.flush() 
    print(f"<<<{total_ways}>>>")


solve_grid_path()