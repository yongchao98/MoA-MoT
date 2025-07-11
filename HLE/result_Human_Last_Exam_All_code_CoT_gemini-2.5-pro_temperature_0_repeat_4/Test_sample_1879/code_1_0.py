def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_TARGET = 4
    U_TARGET = 8
    CONSECUTIVE_MAX = 3

    # DP state: dp[r][u][dir][consecutive]
    # r: number of right moves
    # u: number of up moves
    # dir: 0 for Up, 1 for Right
    # consecutive: number of consecutive moves in 'dir' (1 to CONSECUTIVE_MAX)
    # We use CONSECUTIVE_MAX + 1 for 1-based indexing of consecutive moves.
    dp = [[[[0] * (CONSECUTIVE_MAX + 1) for _ in range(2)] for _ in range(U_TARGET + 1)] for _ in range(R_TARGET + 1)]

    # Base cases: The first step from the origin (0,0)
    # A single step Right to (1,0)
    if R_TARGET > 0:
        dp[1][0][1][1] = 1
    # A single step Up to (0,1)
    if U_TARGET > 0:
        dp[0][1][0][1] = 1

    # Fill the DP table by iterating through each point on the grid
    for r in range(R_TARGET + 1):
        for u in range(U_TARGET + 1):
            # Skip the origin as it's the starting point
            if r == 0 and u == 0:
                continue

            # Calculate ways to reach (r, u) from previous positions

            # Case 1: The last move was RIGHT, coming from (r-1, u)
            if r > 0:
                # To end with exactly 1 Right move, the previous move must have been Up.
                # Sum up all ways to reach (r-1, u) that ended in an Up move.
                ways_from_up = 0
                for k in range(1, CONSECUTIVE_MAX + 1):
                    ways_from_up += dp[r-1][u][0][k]
                dp[r][u][1][1] = ways_from_up

                # To end with k > 1 Right moves, the path to (r-1, u) must have ended with k-1 Right moves.
                for k in range(2, CONSECUTIVE_MAX + 1):
                    dp[r][u][1][k] = dp[r-1][u][1][k-1]

            # Case 2: The last move was UP, coming from (r, u-1)
            if u > 0:
                # To end with exactly 1 Up move, the previous move must have been Right.
                # Sum up all ways to reach (r, u-1) that ended in a Right move.
                ways_from_right = 0
                for k in range(1, CONSECUTIVE_MAX + 1):
                    ways_from_right += dp[r][u-1][1][k]
                dp[r][u][0][1] = ways_from_right

                # To end with k > 1 Up moves, the path to (r, u-1) must have ended with k-1 Up moves.
                for k in range(2, CONSECUTIVE_MAX + 1):
                    dp[r][u][0][k] = dp[r][u-1][0][k-1]

    # The total number of ways is the sum of all valid states at the target (R_TARGET, U_TARGET)
    total_ways = 0
    equation_parts = []
    
    print(f"To find the total number of unique ways, we sum the ways to reach ({R_TARGET},{U_TARGET}) ending with:")
    
    # Sum of ways ending in 1, 2, or 3 'Up' moves
    print("\n- 1, 2, or 3 consecutive 'Up' moves:")
    for k in range(1, CONSECUTIVE_MAX + 1):
        val = dp[R_TARGET][U_TARGET][0][k]
        if val > 0:
            print(f"  Ways ending with {k} 'Up' move(s): {val}")
            total_ways += val
            equation_parts.append(str(val))
            
    # Sum of ways ending in 1, 2, or 3 'Right' moves
    print("\n- 1, 2, or 3 consecutive 'Right' moves:")
    for k in range(1, CONSECUTIVE_MAX + 1):
        val = dp[R_TARGET][U_TARGET][1][k]
        if val > 0:
            print(f"  Ways ending with {k} 'Right' move(s): {val}")
            total_ways += val
            equation_parts.append(str(val))
            
    # Print the final equation and the total
    equation = " + ".join(equation_parts)
    print(f"\nFinal Calculation: {equation} = {total_ways}")

solve_grid_path()
<<<255>>>