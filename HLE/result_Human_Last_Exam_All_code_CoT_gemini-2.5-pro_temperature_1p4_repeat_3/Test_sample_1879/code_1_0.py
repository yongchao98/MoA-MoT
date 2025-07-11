def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with constraints on consecutive moves using dynamic programming.
    """
    # Define the problem parameters
    TARGET_R = 4  # Corresponds to x=4
    TARGET_U = 8  # Corresponds to y=8
    MAX_CONSECUTIVE = 3

    # DP table dimensions: (R+1) x (U+1) x (MAX_CONSECUTIVE+1)
    # ways_r[i][j][k]: number of ways to reach (i, j) with the last k moves being 'Right'
    ways_r = [[[0] * (MAX_CONSECUTIVE + 1) for _ in range(TARGET_U + 1)] for _ in range(TARGET_R + 1)]
    # ways_u[i][j][k]: number of ways to reach (i, j) with the last k moves being 'Up'
    ways_u = [[[0] * (MAX_CONSECUTIVE + 1) for _ in range(TARGET_U + 1)] for _ in range(TARGET_R + 1)]

    # Iterate through the grid, i for right moves, j for up moves
    for i in range(TARGET_R + 1):
        for j in range(TARGET_U + 1):
            if i == 0 and j == 0:
                continue

            # --- Calculate ways to (i, j) ending in Right moves ---
            if i > 0:
                # To end with a single Right (R), the previous path to (i-1, j) must have ended in an Up move.
                total_from_u = sum(ways_u[i-1][j])
                # Base case: The first 'R' move from (0,0) to (1,0). There is 1 way.
                if i == 1 and j == 0:
                    total_from_u = 1
                ways_r[i][j][1] = total_from_u

                # To end with k > 1 Right moves, the path to (i-1, j) must have ended in k-1 Right moves.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    ways_r[i][j][k] = ways_r[i-1][j][k-1]

            # --- Calculate ways to (i, j) ending in Up moves ---
            if j > 0:
                # To end with a single Up (U), the previous path to (i, j-1) must have ended in a Right move.
                total_from_r = sum(ways_r[i][j-1])
                # Base case: The first 'U' move from (0,0) to (0,1). There is 1 way.
                if i == 0 and j == 1:
                    total_from_r = 1
                ways_u[i][j][1] = total_from_r

                # To end with k > 1 Up moves, the path to (i, j-1) must have ended in k-1 Up moves.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    ways_u[i][j][k] = ways_u[i][j-1][k-1]

    # Calculate final results at the target point B(4, 8)
    final_r = [ways_r[TARGET_R][TARGET_U][k] for k in range(1, MAX_CONSECUTIVE + 1)]
    final_u = [ways_u[TARGET_R][TARGET_U][k] for k in range(1, MAX_CONSECUTIVE + 1)]

    total_r = sum(final_r)
    total_u = sum(final_u)
    total_ways = total_r + total_u

    # Print the explanation and final equation
    print(f"Finding unique paths from A(0,0) to B({TARGET_R},{TARGET_U}) with no more than {MAX_CONSECUTIVE} consecutive moves in one direction.\n")
    print(f"Number of paths to ({TARGET_R}, {TARGET_U}) that end with Right (R) moves:")
    for k in range(1, MAX_CONSECUTIVE + 1):
        print(f"  - Ending with {k} consecutive R's: {final_r[k-1]}")

    print(f"\nNumber of paths to ({TARGET_R}, {TARGET_U}) that end with Up (U) moves:")
    for k in range(1, MAX_CONSECUTIVE + 1):
        print(f"  - Ending with {k} consecutive U's: {final_u[k-1]}")

    final_r_str = " + ".join(map(str, final_r))
    final_u_str = " + ".join(map(str, final_u))
    print("\nFinal Calculation:")
    print(f"Total Ways = (Ways ending in R) + (Ways ending in U)")
    print(f"Total Ways = ({final_r_str}) + ({final_u_str})")
    print(f"Total Ways = {total_r} + {total_u}")
    print(f"Total Ways = {total_ways}")
    
    # Final answer in the specified format
    print(f"\n<<<{total_ways}>>>")

solve_grid_path()