def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) with a restriction
    on consecutive moves.
    """
    R_max = 4
    U_max = 8
    K_max = 3  # Maximum number of consecutive steps in the same direction

    # dp[r][u][dir]
    # r: right steps, u: up steps
    # dir: 0 for Right, 1 for Up
    dp = [[[0, 0] for _ in range(U_max + 1)] for _ in range(R_max + 1)]

    # Base case: A virtual starting point at (0,0) to kickstart the DP.
    # This represents one way to be at the start, ready to move Right or Up.
    dp[0][0][0] = 1
    dp[0][0][1] = 1

    # Iterate through the grid to fill the DP table
    for r in range(R_max + 1):
        for u in range(U_max + 1):
            if r == 0 and u == 0:
                continue

            # Calculate ways to reach (r, u) ending with a Right move
            sum_from_u = 0
            for k in range(1, K_max + 1):
                if r - k >= 0:
                    sum_from_u += dp[r - k][u][1]
            dp[r][u][0] = sum_from_u

            # Calculate ways to reach (r, u) ending with an Up move
            sum_from_r = 0
            for k in range(1, K_max + 1):
                if u - k >= 0:
                    sum_from_r += dp[r][u - k][0]
            dp[r][u][1] = sum_from_r

    # The final answer is the sum of ways to reach (R_max, U_max)
    # ending in either a Right or an Up move.
    ways_ending_in_R = dp[R_max][U_max][0]
    ways_ending_in_U = dp[R_max][U_max][1]
    total_ways = ways_ending_in_R + ways_ending_in_U

    print(f"Number of paths ending with a Right move: {ways_ending_in_R}")
    print(f"Number of paths ending with an Up move: {ways_ending_in_U}")
    print(f"Total unique ways: {ways_ending_in_R} + {ways_ending_in_U} = {total_ways}")

solve_grid_path()