def solve_grid_paths():
    """
    Calculates the number of unique ways to move from A(0,0) to B(4,8)
    with the constraint that you cannot move four or more consecutive steps
    in the same direction.
    """
    R_max = 4
    U_max = 8
    consecutive_limit = 3  # Max allowed consecutive moves is 3

    # dp[r][u][dir][k] where:
    # r: number of right moves (0 to 4)
    # u: number of up moves (0 to 8)
    # dir: direction (0 for Right, 1 for Up)
    # k: consecutive moves in that direction, minus 1.
    #    k=0 means 1 move, k=1 means 2 moves, k=2 means 3 moves.
    dp = [[[[0] * consecutive_limit for _ in range(2)] for _ in range(U_max + 1)] for _ in range(R_max + 1)]

    # Iterate through each cell of the grid to calculate paths to it.
    for r in range(R_max + 1):
        for u in range(U_max + 1):
            if r == 0 and u == 0:
                continue

            # --- Calculate ways to arrive at (r, u) with a RIGHT move ---
            if r > 0:
                # Case 1: Arriving with exactly 1 Right move (e.g., ...UR)
                if r == 1 and u == 0:
                    # Path from origin is just 'R', a run of 1.
                    num_ways = 1
                else:
                    # Must come from a path ending in any number of U's at (r-1, u).
                    num_ways = sum(dp[r - 1][u][1])
                dp[r][u][0][0] = num_ways

                # Case 2 & 3: Arriving with 2 or 3 consecutive Right moves.
                for k in range(1, consecutive_limit):
                    # A run of (k+1) R's must follow a run of k R's from the previous cell.
                    dp[r][u][0][k] = dp[r - 1][u][0][k - 1]

            # --- Calculate ways to arrive at (r, u) with an UP move ---
            if u > 0:
                # Case 1: Arriving with exactly 1 Up move (e.g., ...RU)
                if u == 1 and r == 0:
                    # Path from origin is just 'U', a run of 1.
                    num_ways = 1
                else:
                    # Must come from a path ending in any number of R's at (r, u-1).
                    num_ways = sum(dp[r][u - 1][0])
                dp[r][u][1][0] = num_ways
                
                # Case 2 & 3: Arriving with 2 or 3 consecutive Up moves.
                for k in range(1, consecutive_limit):
                    # A run of (k+1) U's must follow a run of k U's from the previous cell.
                    dp[r][u][1][k] = dp[r][u - 1][1][k - 1]

    # Collect the final results for the destination B(4,8).
    ways_ending_R1 = dp[R_max][U_max][0][0]
    ways_ending_R2 = dp[R_max][U_max][0][1]
    ways_ending_R3 = dp[R_max][U_max][0][2]
    
    ways_ending_U1 = dp[R_max][U_max][1][0]
    ways_ending_U2 = dp[R_max][U_max][1][1]
    ways_ending_U3 = dp[R_max][U_max][1][2]

    total_ways = (ways_ending_R1 + ways_ending_R2 + ways_ending_R3 +
                  ways_ending_U1 + ways_ending_U2 + ways_ending_U3)

    # Print the breakdown of the final count as an equation.
    print("The final answer is the sum of paths ending in different states:")
    print(f"{ways_ending_R1} (ends in R) + {ways_ending_R2} (ends in RR) + {ways_ending_R3} (ends in RRR) + "
          f"{ways_ending_U1} (ends in U) + {ways_ending_U2} (ends in UU) + {ways_ending_U3} (ends in UUU) "
          f"= {total_ways}")

solve_grid_paths()
<<<254>>>