def solve_grid_walk():
    """
    Calculates the number of unique ways to walk from (0,0) to (4,8)
    with the given constraints using dynamic programming.
    """
    R_MAX = 4
    U_MAX = 8
    # The maximum number of consecutive steps allowed is 3.
    K_MAX = 3

    # dp[r][u][dir][k]: ways to get to (r, u)
    # r: number of right moves
    # u: number of up moves
    # dir: 0 for Right, 1 for Up
    # k: number of consecutive moves in `dir` (1 to K_MAX)
    # We use k+1 for array size to allow 1-based indexing for k.
    dp = [[[[0] * (K_MAX + 1) for _ in range(2)] for _ in range(U_MAX + 1)] for _ in range(R_MAX + 1)]

    # Base cases: From (0,0), we can only go to (1,0) or (0,1).
    # One way to reach (1,0) ending with 1 Right move.
    if R_MAX > 0:
        dp[1][0][0][1] = 1
    # One way to reach (0,1) ending with 1 Up move.
    if U_MAX > 0:
        dp[0][1][1][1] = 1

    # Fill the DP table by iterating through the grid
    for r in range(R_MAX + 1):
        for u in range(U_MAX + 1):
            # Skip the origin
            if r == 0 and u == 0:
                continue

            # Calculate ways to reach (r, u)
            
            # Case 1: Last move was Right (came from (r-1, u))
            if r > 0:
                # Option A: Previous move was also Right.
                # Add a Right move to a sequence of k-1 Right moves.
                # This is only possible if k-1 < K_MAX.
                for k in range(2, K_MAX + 1):
                    dp[r][u][0][k] += dp[r-1][u][0][k-1]
                
                # Option B: Previous move was Up.
                # Start a new sequence of Right moves (k=1).
                # This can come from any valid sequence of Up moves ending at (r-1, u).
                for k in range(1, K_MAX + 1):
                    dp[r][u][0][1] += dp[r-1][u][1][k]

            # Case 2: Last move was Up (came from (r, u-1))
            if u > 0:
                # Option A: Previous move was also Up.
                # Add an Up move to a sequence of k-1 Up moves.
                for k in range(2, K_MAX + 1):
                    dp[r][u][1][k] += dp[r][u-1][1][k-1]

                # Option B: Previous move was Right.
                # Start a new sequence of Up moves (k=1).
                for k in range(1, K_MAX + 1):
                    dp[r][u][1][1] += dp[r][u-1][0][k]

    # The final answer is the sum of all ways to reach the destination (R_MAX, U_MAX)
    # ending in any valid sequence of moves.
    final_counts = []
    
    print("Number of ways to reach B(4,8):")
    # Sum ways ending in Right moves
    for k in range(1, K_MAX + 1):
        count = dp[R_MAX][U_MAX][0][k]
        final_counts.append(count)
        print(f"Ending with {k} consecutive Right moves: {count}")

    # Sum ways ending in Up moves
    for k in range(1, K_MAX + 1):
        count = dp[R_MAX][U_MAX][1][k]
        final_counts.append(count)
        print(f"Ending with {k} consecutive Up moves: {count}")
    
    total_ways = sum(final_counts)
    
    equation = " + ".join(map(str, final_counts))
    print(f"\nTotal unique ways = {equation} = {total_ways}")
    
    return total_ways

# Execute the function and get the final answer
final_answer = solve_grid_walk()
print(f"<<<{final_answer}>>>")