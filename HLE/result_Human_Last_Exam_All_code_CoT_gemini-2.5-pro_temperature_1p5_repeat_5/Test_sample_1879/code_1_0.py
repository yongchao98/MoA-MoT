def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid,
    with the constraint that no more than 3 consecutive steps in the same direction
    are allowed.
    """
    max_R = 4
    max_U = 8
    # The maximum number of allowed consecutive moves is 3.
    # We will use an array of size 3, where index k corresponds to k+1 moves.
    # index 0 -> 1 move, index 1 -> 2 moves, index 2 -> 3 moves.
    max_consecutive_len = 3

    # R_counts[r][u][k] = # of paths to (r,u) ending in k+1 'R's
    R_counts = [[[0] * max_consecutive_len for _ in range(max_U + 1)] for _ in range(max_R + 1)]
    # U_counts[r][u][k] = # of paths to (r,u) ending in k+1 'U's
    U_counts = [[[0] * max_consecutive_len for _ in range(max_U + 1)] for _ in range(max_R + 1)]

    # Initialize the process. After one step, we are at (1,0) or (0,1).
    # There is 1 way to get to (1,0) ending in one 'R'.
    if max_R >= 1:
        R_counts[1][0][0] = 1
    # There is 1 way to get to (0,1) ending in one 'U'.
    if max_U >= 1:
        U_counts[0][1][0] = 1

    # Iterate through the grid and fill the DP tables.
    for r in range(max_R + 1):
        for u in range(max_U + 1):
            # Skip the points we already handled as base cases
            if (r == 0 and u == 0) or (r == 1 and u == 0) or (r == 0 and u == 1):
                continue

            # Calculate paths to (r,u) ending in 'R's.
            # This requires coming from (r-1, u).
            if r > 0:
                # To end in exactly one 'R', the previous path must have ended in 'U'.
                total_from_U = sum(U_counts[r-1][u])
                R_counts[r][u][0] = total_from_U

                # To end in k+1 'R's, the previous path must have ended in k 'R's.
                for k in range(1, max_consecutive_len):
                    R_counts[r][u][k] = R_counts[r-1][u][k-1]

            # Calculate paths to (r,u) ending in 'U's.
            # This requires coming from (r, u-1).
            if u > 0:
                # To end in exactly one 'U', the previous path must have ended in 'R'.
                total_from_R = sum(R_counts[r][u-1])
                U_counts[r][u][0] = total_from_R

                # To end in k+1 'U's, the previous path must have ended in k 'U's.
                for k in range(1, max_consecutive_len):
                    U_counts[r][u][k] = U_counts[r][u-1][k-1]

    # The final answer is the sum of all paths to the destination (4,8).
    final_R_ways = sum(R_counts[max_R][max_U])
    final_U_ways = sum(U_counts[max_R][max_U])
    total_ways = final_R_ways + final_U_ways

    print(f"The calculation is based on the final move to reach ({max_R},{max_U}):")
    print(f"Number of ways ending with a 'Right' move: {final_R_ways}")
    print(f"Number of ways ending with an 'Up' move: {final_U_ways}")
    print("\nThe final equation for the total number of ways is:")
    print(f"{final_R_ways} + {final_U_ways} = {total_ways}")
    print(f"\nThus, there are {total_ways} unique ways to move from A(0,0) to B(4,8).")
    
    return total_ways

# Run the solver and store the answer
result = solve_grid_path()
print(f'<<<{result}>>>')
