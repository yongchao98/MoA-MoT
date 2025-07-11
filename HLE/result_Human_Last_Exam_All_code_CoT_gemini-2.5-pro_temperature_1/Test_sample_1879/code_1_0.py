def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_DEST = 4
    U_DEST = 8
    MAX_CONSECUTIVE = 3

    # DP table: dp[r][u][dir][k]
    # dir: 0 for Right, 1 for Up
    # k: number of consecutive steps (indices 1, 2, 3)
    dp = [[[[0] * (MAX_CONSECUTIVE + 1) for _ in range(2)] for _ in range(U_DEST + 1)] for _ in range(R_DEST + 1)]

    # Base cases: A single step from the origin (0,0)
    # One step Right to (1,0)
    if R_DEST > 0:
        dp[1][0][0][1] = 1
    # One step Up to (0,1)
    if U_DEST > 0:
        dp[0][1][1][1] = 1

    # Fill the DP table by iterating through the grid
    for r in range(R_DEST + 1):
        for u in range(U_DEST + 1):
            # Skip the origin and the base cases already set
            if r == 0 and u == 0: continue
            if r == 1 and u == 0: continue
            if r == 0 and u == 1: continue

            # --- Calculate ways to reach (r, u) ---

            # Case 1: The last move was RIGHT (coming from r-1, u)
            if r > 0:
                # To end with k=1 RIGHT move, the previous move must be UP.
                # Sum of all ways to reach (r-1, u) ending in UP.
                ways_ending_in_up = dp[r-1][u][1][1] + dp[r-1][u][1][2] + dp[r-1][u][1][3]
                dp[r][u][0][1] = ways_ending_in_up

                # To end with k>1 RIGHT moves, extend a path ending in k-1 RIGHT moves.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    dp[r][u][0][k] = dp[r-1][u][0][k-1]
            
            # Case 2: The last move was UP (coming from r, u-1)
            if u > 0:
                # To end with k=1 UP move, the previous move must be RIGHT.
                # Sum of all ways to reach (r, u-1) ending in RIGHT.
                ways_ending_in_right = dp[r][u-1][0][1] + dp[r][u-1][0][2] + dp[r][u-1][0][3]
                dp[r][u][1][1] = ways_ending_in_right
                
                # To end with k>1 UP moves, extend a path ending in k-1 UP moves.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    dp[r][u][1][k] = dp[r][u-1][1][k-1]

    # The final answer is the sum of all ways to reach the destination (R_DEST, U_DEST)
    final_dp_right = dp[R_DEST][U_DEST][0]
    final_dp_up = dp[R_DEST][U_DEST][1]

    total_ways = sum(final_dp_right) + sum(final_dp_up)
    
    print("Ways to reach (4,8) ending with:")
    print(f"- 1 Right step: {final_dp_right[1]}")
    print(f"- 2 Right steps: {final_dp_right[2]}")
    print(f"- 3 Right steps: {final_dp_right[3]}")
    print(f"- 1 Up step:    {final_dp_up[1]}")
    print(f"- 2 Up steps:    {final_dp_up[2]}")
    print(f"- 3 Up steps:    {final_dp_up[3]}")
    print("\nFinal calculation:")
    
    # Print the equation as requested
    numbers_to_sum = [
        final_dp_right[1], final_dp_right[2], final_dp_right[3],
        final_dp_up[1], final_dp_up[2], final_dp_up[3]
    ]
    equation = " + ".join(map(str, numbers_to_sum))
    print(f"{equation} = {total_ways}")

    print(f"\nTotal unique ways: {total_ways}")
    print(f"<<<{total_ways}>>>")

solve_grid_path()