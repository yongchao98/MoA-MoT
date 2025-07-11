def solve_grid_path():
    """
    Calculates the number of unique ways to walk from A(0,0) to B(4,8)
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_TARGET = 4
    U_TARGET = 8
    MAX_CONSECUTIVE = 3

    # dp[r][u][dir][k]
    # r: right steps, u: up steps
    # dir: 0 for Right, 1 for Up
    # k: number of consecutive steps - 1 (k=0 for 1, k=1 for 2, k=2 for 3)
    dp = [[[[0] * MAX_CONSECUTIVE for _ in range(2)] for _ in range(U_TARGET + 1)] for _ in range(R_TARGET + 1)]

    # We iterate by the total number of steps to ensure dependencies are met.
    for s in range(1, R_TARGET + U_TARGET + 1):
        for r in range(s + 1):
            u = s - r
            if r > R_TARGET or u > U_TARGET:
                continue

            # --- Calculate ways ending in Right moves to reach (r, u) ---
            if r > 0:
                # Case 1: Ending with a single 'R'. This move must follow a 'U' move,
                # or it's the very first step from the origin.
                if r == 1 and u == 0:
                    ways_from_u = 1  # The path is just "R" from the origin.
                else:
                    # Sum of all ways to get to (r-1, u) that ended in a 'U' move.
                    ways_from_u = sum(dp[r - 1][u][1])
                dp[r][u][0][0] = ways_from_u

                # Case 2: Ending with more than one 'R'.
                for k in range(1, MAX_CONSECUTIVE):
                    # A path ending in k+1 'R's comes from a path ending in k 'R's at (r-1, u).
                    dp[r][u][0][k] = dp[r - 1][u][0][k - 1]

            # --- Calculate ways ending in Up moves to reach (r, u) ---
            if u > 0:
                # Case 1: Ending with a single 'U'.
                if r == 0 and u == 1:
                    ways_from_r = 1  # The path is just "U" from the origin.
                else:
                    # Sum of all ways to get to (r, u-1) that ended in an 'R' move.
                    ways_from_r = sum(dp[r][u - 1][0])
                dp[r][u][1][0] = ways_from_r

                # Case 2: Ending with more than one 'U'.
                for k in range(1, MAX_CONSECUTIVE):
                    # A path ending in k+1 'U's comes from a path ending in k 'U's at (r, u-1).
                    dp[r][u][1][k] = dp[r][u - 1][1][k - 1]

    # The final answer is the sum of all possible ways to arrive at the destination.
    final_r_ways = dp[R_TARGET][U_TARGET][0]
    final_u_ways = dp[R_TARGET][U_TARGET][1]

    total_ways = sum(final_r_ways) + sum(final_u_ways)
    
    print("Number of paths ending with:")
    print(f"- 1 Right move: {final_r_ways[0]}")
    print(f"- 2 Right moves: {final_r_ways[1]}")
    print(f"- 3 Right moves: {final_r_ways[2]}")
    print(f"- 1 Up move: {final_u_ways[0]}")
    print(f"- 2 Up moves: {final_u_ways[1]}")
    print(f"- 3 Up moves: {final_u_ways[2]}")
    
    equation_parts = final_r_ways + final_u_ways
    equation_str = " + ".join(map(str, equation_parts))
    
    print(f"\nTotal unique ways = {equation_str} = {total_ways}")
    return total_ways

result = solve_grid_path()
print(f"<<<{result}>>>")