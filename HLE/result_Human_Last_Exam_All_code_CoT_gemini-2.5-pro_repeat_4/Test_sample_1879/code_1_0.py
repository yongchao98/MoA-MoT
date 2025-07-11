def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with specific movement constraints using dynamic programming.
    """
    target_r = 4
    target_u = 8
    # The maximum number of consecutive steps allowed is 3.
    max_consecutive = 3

    # DP table: dp[r][u][dir][k]
    # r: 0 to target_r (number of right moves)
    # u: 0 to target_u (number of up moves)
    # dir: 0 for Right, 1 for Up
    # k: 0 to max_consecutive-1 (representing 1 to max_consecutive steps)
    dp = [[[[0] * max_consecutive for _ in range(2)] for _ in range(target_u + 1)] for _ in range(target_r + 1)]

    # Base cases: Starting from (0,0)
    # One 'Right' move to (1,0)
    if 1 <= target_r:
        dp[1][0][0][0] = 1
    # One 'Up' move to (0,1)
    if 1 <= target_u:
        dp[0][1][1][0] = 1

    # Fill the DP table by iterating through the grid
    for r in range(target_r + 1):
        for u in range(target_u + 1):
            # Paths ending at (r, u) with a RIGHT move (must come from (r-1, u))
            if r > 0:
                # To end with 1 R move, the previous move must have been U
                total_from_up = 0
                for k in range(max_consecutive):
                    total_from_up += dp[r - 1][u][1][k]
                dp[r][u][0][0] = total_from_up

                # To end with k+1 R moves, the path to (r-1, u) must have ended in k R moves
                for k in range(1, max_consecutive):
                    dp[r][u][0][k] = dp[r - 1][u][0][k - 1]

            # Paths ending at (r, u) with an UP move (must come from (r, u-1))
            if u > 0:
                # To end with 1 U move, the previous move must have been R
                total_from_right = 0
                for k in range(max_consecutive):
                    total_from_right += dp[r][u - 1][0][k]
                dp[r][u][1][0] = total_from_right

                # To end with k+1 U moves, the path to (r, u-1) must have ended in k U moves
                for k in range(1, max_consecutive):
                    dp[r][u][1][k] = dp[r][u - 1][1][k - 1]

    # The final answer is the sum of all ways to reach the target (4,8)
    total_ways = 0
    components = []
    
    # Sum paths ending in Right moves
    for k in range(max_consecutive):
        ways = dp[target_r][target_u][0][k]
        print(f"Paths to ({target_r}, {target_u}) ending with {k+1} consecutive Right moves: {ways}")
        total_ways += ways
        components.append(str(ways))

    # Sum paths ending in Up moves
    for k in range(max_consecutive):
        ways = dp[target_r][target_u][1][k]
        print(f"Paths to ({target_r}, {target_u}) ending with {k+1} consecutive Up moves: {ways}")
        total_ways += ways
        components.append(str(ways))
        
    # Print the final equation
    print("\nFinal Calculation:")
    equation = " + ".join(components)
    print(f"{equation} = {total_ways}")
    print(f"\nTotal unique ways to move from A(0,0) to B({target_r},{target_u}) are: {total_ways}")

solve_grid_path()
<<<270>>>