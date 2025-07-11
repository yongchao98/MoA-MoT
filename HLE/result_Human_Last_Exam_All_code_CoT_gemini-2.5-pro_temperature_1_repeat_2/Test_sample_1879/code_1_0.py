def solve_grid_walk():
    """
    Calculates the number of unique paths from (0,0) to (R,U) with a constraint
    on the maximum number of consecutive moves in the same direction.
    """
    R, U = 4, 8
    # The maximum number of consecutive steps allowed is 3.
    MAX_CONSECUTIVE = 3

    # dp_h[r][u][k] = ways to get to (r,u) ending in k consecutive R's.
    # dp_v[r][u][k] = ways to get to (r,u) ending in k consecutive U's.
    # k is 1-based, so we use size MAX_CONSECUTIVE + 1.
    # Indices: r for right moves, u for up moves, k for consecutive count.
    dp_h = [[[0] * (MAX_CONSECUTIVE + 1) for _ in range(U + 1)] for _ in range(R + 1)]
    dp_v = [[[0] * (MAX_CONSECUTIVE + 1) for _ in range(U + 1)] for _ in range(R + 1)]

    # Base case: At the start (0,0), we can make one move.
    # A single 'R' move to (1,0) ends in 1 consecutive 'R'.
    if R > 0:
        dp_h[1][0][1] = 1
    # A single 'U' move to (0,1) ends in 1 consecutive 'U'.
    if U > 0:
        dp_v[0][1][1] = 1

    # Fill the DP tables by iterating through total steps s = r + u
    for s in range(2, R + U + 1):
        for r in range(s + 1):
            u = s - r
            if r > R or u > U:
                continue

            # Calculate ways to reach (r, u) ending in a Right move
            if r > 0:
                # To end with 1 'R', the previous move must be 'U'.
                # Sum up all ways to get to (r-1, u) that ended in 'U'.
                total_v_prev = sum(dp_v[r - 1][u])
                dp_h[r][u][1] = total_v_prev
                
                # To end with k 'R's (k>1), the previous state must have ended in k-1 'R's.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    dp_h[r][u][k] = dp_h[r - 1][u][k - 1]

            # Calculate ways to reach (r, u) ending in an Up move
            if u > 0:
                # To end with 1 'U', the previous move must be 'R'.
                # Sum up all ways to get to (r, u-1) that ended in 'R'.
                total_h_prev = sum(dp_h[r][u - 1])
                dp_v[r][u][1] = total_h_prev

                # To end with k 'U's (k>1), the previous state must have ended in k-1 'U's.
                for k in range(2, MAX_CONSECUTIVE + 1):
                    dp_v[r][u][k] = dp_v[r][u - 1][k - 1]

    # The final answer is the sum of all ways to reach (R, U).
    final_h = dp_h[R][U]
    final_v = dp_v[R][U]
    total_ways = sum(final_h) + sum(final_v)
    
    print(f"To reach ({R},{U}), there are:")
    print(f" {sum(final_h)} ways ending in a Right move:")
    print(f"  Ending in 1 R: {final_h[1]}")
    print(f"  Ending in 2 R: {final_h[2]}")
    print(f"  Ending in 3 R: {final_h[3]}")
    print(f" {sum(final_v)} ways ending in an Up move:")
    print(f"  Ending in 1 U: {final_v[1]}")
    print(f"  Ending in 2 U: {final_v[2]}")
    print(f"  Ending in 3 U: {final_v[3]}")
    print("-" * 20)
    print("The final equation for the total number of ways is:")
    print(f"{final_h[1]} + {final_h[2]} + {final_h[3]} + {final_v[1]} + {final_v[2]} + {final_v[3]} = {total_ways}")

solve_grid_walk()