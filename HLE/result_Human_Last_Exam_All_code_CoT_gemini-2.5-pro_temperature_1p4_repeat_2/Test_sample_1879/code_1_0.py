def solve_grid_path():
    """
    Calculates the number of unique paths on a 2D grid from (0,0) to (4,8)
    with specific movement constraints using dynamic programming.
    """
    R_max = 4
    U_max = 8

    # DP tables: R_paths[r][u][k] is the number of ways to reach (r,u)
    # ending with exactly k 'Right' moves. k can be 1, 2, or 3.
    # We use size k=4 (indices 0-3) and leave index 0 unused for clarity.
    R_paths = [[[0] * 4 for _ in range(U_max + 1)] for _ in range(R_max + 1)]
    U_paths = [[[0] * 4 for _ in range(U_max + 1)] for _ in range(R_max + 1)]

    # Base cases for the first step from (0,0)
    # Path 'R' to (1,0)
    if R_max > 0:
        R_paths[1][0][1] = 1
    # Path 'U' to (0,1)
    if U_max > 0:
        U_paths[0][1][1] = 1

    # Fill the DP tables by iterating through the grid
    for r in range(R_max + 1):
        for u in range(U_max + 1):
            # Skip the origin, it has no paths leading to it
            if r == 0 and u == 0:
                continue

            # Calculate paths to (r, u)
            
            # Paths ending in 'Right' moves
            if r > 0:
                # Paths ending in exactly one 'R' must come from a path ending in 'U' at (r-1, u)
                sum_u_prev = sum(U_paths[r - 1][u])
                R_paths[r][u][1] = sum_u_prev
                # Paths ending in 2 or 3 'R's must come from a path ending in 1 or 2 'R's at (r-1, u)
                R_paths[r][u][2] = R_paths[r - 1][u][1]
                R_paths[r][u][3] = R_paths[r - 1][u][2]

            # Paths ending in 'Up' moves
            if u > 0:
                # Paths ending in exactly one 'U' must come from a path ending in 'R' at (r, u-1)
                sum_r_prev = sum(R_paths[r][u - 1])
                U_paths[r][u][1] = sum_r_prev
                # Paths ending in 2 or 3 'U's must come from a path ending in 1 or 2 'U's at (r, u-1)
                U_paths[r][u][2] = U_paths[r][u - 1][1]
                U_paths[r][u][3] = U_paths[r][u - 1][2]

    # The final answer is the sum of all paths reaching the destination (R_max, U_max)
    final_R_paths = R_paths[R_max][U_max]
    final_U_paths = U_paths[R_max][U_max]
    
    total_R = sum(final_R_paths)
    total_U = sum(final_U_paths)
    total_ways = total_R + total_U

    print(f"Paths to ({R_max}, {U_max}) ending with Right moves:")
    print(f" - Ending in 1 'R': {final_R_paths[1]}")
    print(f" - Ending in 2 'R's: {final_R_paths[2]}")
    print(f" - Ending in 3 'R's: {final_R_paths[3]}")
    print(f"Total paths ending in R = {final_R_paths[1]} + {final_R_paths[2]} + {final_R_paths[3]} = {total_R}\n")

    print(f"Paths to ({R_max}, {U_max}) ending with Up moves:")
    print(f" - Ending in 1 'U': {final_U_paths[1]}")
    print(f" - Ending in 2 'U's: {final_U_paths[2]}")
    print(f" - Ending in 3 'U's: {final_U_paths[3]}")
    print(f"Total paths ending in U = {final_U_paths[1]} + {final_U_paths[2]} + {final_U_paths[3]} = {total_U}\n")

    print(f"Total unique ways = {total_R} + {total_U} = {total_ways}")
    
    return total_ways

total_ways = solve_grid_path()
print(f"<<<{total_ways}>>>")