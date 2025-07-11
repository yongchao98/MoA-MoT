def solve_grid_paths():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint that no four or more consecutive steps in the same
    direction are allowed.
    """
    R_max = 4
    U_max = 8
    # K is the maximum number of consecutive moves allowed.
    # The problem states "four or more" are not allowed, so max is 3.
    K = 3

    # Initialize DP table: dp[r][u][dir][len]
    # Dimensions: (R_max+1) x (U_max+1) x 2 (directions) x (K+1) (lengths 0-3)
    # dir: 0 for Right, 1 for Up
    # len: consecutive moves in that direction (1-indexed)
    dp = [[[([0] * (K + 1)) for _ in range(2)] for _ in range(U_max + 1)] for _ in range(R_max + 1)]

    # Loop through the grid to fill the DP table
    for r in range(R_max + 1):
        for u in range(U_max + 1):
            if r == 0 and u == 0:
                continue

            # Calculate ways to get to (r, u)
            
            # Case 1: The last move was Right (came from (r-1, u))
            if r > 0:
                # Subcase 1.1: The path to (r-1, u) ended in Up moves.
                # This starts a new sequence of Right moves (length 1).
                # We need to sum up all ways to get to (r-1, u) that ended in 'Up'.
                if r == 1 and u == 0:
                    # Base case: Path is just 'R' from the origin (0,0).
                    # We consider the path from origin as 1 way.
                    ways_from_up_path = 1
                else:
                    ways_from_up_path = sum(dp[r - 1][u][1][k] for k in range(1, K + 1))
                dp[r][u][0][1] = ways_from_up_path

                # Subcase 1.2: The path to (r-1, u) ended in Right moves.
                # This extends the sequence of Right moves.
                for k in range(2, K + 1):
                    dp[r][u][0][k] = dp[r - 1][u][0][k - 1]

            # Case 2: The last move was Up (came from (r, u-1))
            if u > 0:
                # Subcase 2.1: The path to (r, u-1) ended in Right moves.
                # This starts a new sequence of Up moves (length 1).
                if r == 0 and u == 1:
                    # Base case: Path is just 'U' from the origin (0,0).
                    ways_from_right_path = 1
                else:
                    ways_from_right_path = sum(dp[r][u - 1][0][k] for k in range(1, K + 1))
                dp[r][u][1][1] = ways_from_right_path
                
                # Subcase 2.2: The path to (r, u-1) ended in Up moves.
                # This extends the sequence of Up moves.
                for k in range(2, K + 1):
                    dp[r][u][1][k] = dp[r][u - 1][1][k - 1]

    # The final answer is the sum of all ways to reach the destination (R_max, U_max)
    # ending in any allowed state.
    
    ways_ending_R = [dp[R_max][U_max][0][k] for k in range(1, K + 1)]
    ways_ending_U = [dp[R_max][U_max][1][k] for k in range(1, K + 1)]

    total_ways = sum(ways_ending_R) + sum(ways_ending_U)
    
    print("Paths to (4,8) can end with:")
    print(f"- 1 Right move: {ways_ending_R[0]} ways")
    print(f"- 2 consecutive Right moves: {ways_ending_R[1]} ways")
    print(f"- 3 consecutive Right moves: {ways_ending_R[2]} ways")
    print(f"- 1 Up move: {ways_ending_U[0]} ways")
    print(f"- 2 consecutive Up moves: {ways_ending_U[1]} ways")
    print(f"- 3 consecutive Up moves: {ways_ending_U[2]} ways")
    print("\nFinal Calculation:")
    
    r_sum_str = " + ".join(map(str, ways_ending_R))
    u_sum_str = " + ".join(map(str, ways_ending_U))
    
    print(f"Total ways = (Ways ending in R) + (Ways ending in U)")
    print(f"Total ways = ({r_sum_str}) + ({u_sum_str})")
    print(f"Total ways = {sum(ways_ending_R)} + {sum(ways_ending_U)} = {total_ways}")
    
# Execute the function to find the solution.
solve_grid_paths()