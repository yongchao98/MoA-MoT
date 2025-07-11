def solve_grid_path_constraint():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid,
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_target = 4
    U_target = 8
    # The constraint is "four or more", so the max allowed is 3.
    consecutive_limit = 3

    # dp_R[r][u]: number of paths to (r, u) ending in a block of Right moves.
    # dp_U[r][u]: number of paths to (r, u) ending in a block of Up moves.
    dp_R = [[0] * (U_target + 1) for _ in range(R_target + 1)]
    dp_U = [[0] * (U_target + 1) for _ in range(R_target + 1)]

    # We use a "phantom" base case at (0,0) to start the recurrence.
    # This represents an empty path at the origin, from which we can start
    # a sequence of moves with either R or U.
    # dp_U[0][0] = 1 acts as the base for any path starting with R.
    # dp_R[0][0] = 1 acts as the base for any path starting with U.
    dp_U[0][0] = 1
    dp_R[0][0] = 1

    # Iterate through each point on the grid to calculate the number of paths.
    for r in range(R_target + 1):
        for u in range(U_target + 1):
            if r == 0 and u == 0:
                continue

            # Calculate paths to (r, u) ending in a block of Right moves.
            # This is the sum of paths coming from (r-k, u) that ended in Up moves.
            paths_ending_in_R = 0
            for k in range(1, consecutive_limit + 1):
                if r - k >= 0:
                    paths_ending_in_R += dp_U[r - k][u]
            dp_R[r][u] = paths_ending_in_R

            # Calculate paths to (r, u) ending in a block of Up moves.
            # This is the sum of paths coming from (r, u-k) that ended in Right moves.
            paths_ending_in_U = 0
            for k in range(1, consecutive_limit + 1):
                if u - k >= 0:
                    paths_ending_in_U += dp_R[r][u - k]
            dp_U[r][u] = paths_ending_in_U

    # The final answer is the sum of paths to the target B(4,8)
    # that end in either a Right move or an Up move.
    paths_ending_R_at_B = dp_R[R_target][U_target]
    paths_ending_U_at_B = dp_U[R_target][U_target]
    total_ways = paths_ending_R_at_B + paths_ending_U_at_B

    print(f"Number of paths to B({R_target}, {U_target}) that end with a Right move = {paths_ending_R_at_B}")
    print(f"Number of paths to B({R_target}, {U_target}) that end with an Up move = {paths_ending_U_at_B}")
    print(f"Total Unique Ways = {paths_ending_R_at_B} + {paths_ending_U_at_B} = {total_ways}")

solve_grid_path_constraint()
<<<168>>>