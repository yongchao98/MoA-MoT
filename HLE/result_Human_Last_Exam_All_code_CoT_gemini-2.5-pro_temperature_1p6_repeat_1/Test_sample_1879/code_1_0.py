def solve_grid_path():
    """
    Calculates the number of unique ways to move from (0,0) to (4,8) on a 2D grid
    with the constraint of not having 4 or more consecutive steps in the same direction.
    """
    R_MOVES = 4
    U_MOVES = 8

    # dp[i][j] will store a tuple: (paths_ending_with_R, paths_ending_with_U)
    # paths_ending_with_R is a list [r1, r2, r3] where r_k is the number of paths ending with k R's.
    # paths_ending_with_U is a list [u1, u2, u3] where u_k is the number of paths ending with k U's.
    dp = [[([0, 0, 0], [0, 0, 0]) for _ in range(U_MOVES + 1)] for _ in range(R_MOVES + 1)]

    # Base case: At (0,0), we imagine a path of length 0.
    # We can start by moving Right to (1,0) or Up to (0,1).

    # Initialize paths along the axes
    # Paths to (i, 0) can only consist of 'R' moves
    for i in range(1, R_MOVES + 1):
        if i < 4:
            r_counts = [0, 0, 0]
            r_counts[i - 1] = 1
            dp[i][0] = (r_counts, [0, 0, 0])

    # Paths to (0, j) can only consist of 'U' moves
    for j in range(1, U_MOVES + 1):
        if j < 4:
            u_counts = [0, 0, 0]
            u_counts[j - 1] = 1
            dp[0][j] = ([0, 0, 0], u_counts)

    # Fill the DP table
    for i in range(1, R_MOVES + 1):
        for j in range(1, U_MOVES + 1):
            # Calculate paths to (i,j) ending with R moves
            # These must come from (i-1, j)
            prev_u_counts = dp[i-1][j][1]
            prev_r_counts = dp[i-1][j][0]
            r1 = sum(prev_u_counts)
            r2 = prev_r_counts[0]
            r3 = prev_r_counts[1]
            
            # Calculate paths to (i,j) ending with U moves
            # These must come from (i, j-1)
            prev_u_counts_up = dp[i][j-1][1]
            prev_r_counts_up = dp[i][j-1][0]
            u1 = sum(prev_r_counts_up)
            u2 = prev_u_counts_up[0]
            u3 = prev_u_counts_up[1]

            dp[i][j] = ([r1, r2, r3], [u1, u2, u3])
    
    final_r_counts = dp[R_MOVES][U_MOVES][0]
    final_u_counts = dp[R_MOVES][U_MOVES][1]
    
    total_ways = sum(final_r_counts) + sum(final_u_counts)

    # Output the final summation
    r_sum_str = " + ".join(map(str, final_r_counts))
    u_sum_str = " + ".join(map(str, final_u_counts))
    
    print(f"The number of valid paths is the sum of paths ending in R and paths ending in U.")
    print(f"Paths ending in R (1, 2, or 3 consecutive): {r_sum_str}")
    print(f"Paths ending in U (1, 2, or 3 consecutive): {u_sum_str}")
    print(f"Final calculation: {r_sum_str} + {u_sum_str} = {total_ways}")

solve_grid_path()