def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_MAX = 4
    U_MAX = 8
    CONSECUTIVE_MAX = 3

    # Directions: 0 for Right, 1 for Up
    R_DIR, U_DIR = 0, 1

    # DP table: dp[i][j][dir][k]
    # i: x-coord (Right moves), j: y-coord (Up moves)
    # dir: last move direction, k: number of consecutive moves in dir
    dp = [[[[0] * (CONSECUTIVE_MAX + 1) for _ in range(2)] for _ in range(U_MAX + 1)] for _ in range(R_MAX + 1)]

    # Base Cases: Paths along the axes
    # Paths to (i, 0) must be R, RR, or RRR
    for i in range(1, CONSECUTIVE_MAX + 1):
        if i <= R_MAX:
            dp[i][0][R_DIR][i] = 1

    # Paths to (0, j) must be U, UU, or UUU
    for j in range(1, CONSECUTIVE_MAX + 1):
        if j <= U_MAX:
            dp[0][j][U_DIR][j] = 1

    # Fill the DP table using the recurrence relations
    for i in range(R_MAX + 1):
        for j in range(U_MAX + 1):
            # Skip base cases already handled on the axes
            if i == 0 or j == 0:
                continue

            # Calculate ways to reach (i,j) ending with a Right move
            # This requires coming from (i-1, j)
            v_sum_prev = sum(dp[i-1][j][U_DIR][k] for k in range(1, CONSECUTIVE_MAX + 1))
            dp[i][j][R_DIR][1] = v_sum_prev
            for k in range(2, CONSECUTIVE_MAX + 1):
                dp[i][j][R_DIR][k] = dp[i-1][j][R_DIR][k-1]

            # Calculate ways to reach (i,j) ending with an Up move
            # This requires coming from (i, j-1)
            h_sum_prev = sum(dp[i][j-1][R_DIR][k] for k in range(1, CONSECUTIVE_MAX + 1))
            dp[i][j][U_DIR][1] = h_sum_prev
            for k in range(2, CONSECUTIVE_MAX + 1):
                dp[i][j][U_DIR][k] = dp[i][j-1][U_DIR][k-1]

    # The final answer is the sum of all ways to reach the destination (R_MAX, U_MAX)
    ways_ending_r = sum(dp[R_MAX][U_MAX][R_DIR])
    ways_ending_u = sum(dp[R_MAX][U_MAX][U_DIR])
    total_ways = ways_ending_r + ways_ending_u

    print(f"Number of paths to (4,8) ending with a Right move: {ways_ending_r}")
    print(f"Number of paths to (4,8) ending with an Up move: {ways_ending_u}")
    print("\nFinal Equation:")
    print(f"{ways_ending_r} + {ways_ending_u} = {total_ways}")
    print(f"\nTotal unique ways from A(0,0) to B(4,8) are: {total_ways}")
    return total_ways

# Run the solver
final_answer = solve_grid_path()
# The final answer format is specified in the problem description.
# The code above prints the required components and the final sum.
# The numeric answer itself is returned by the function.
# e.g., print(f"<<<{final_answer}>>>")