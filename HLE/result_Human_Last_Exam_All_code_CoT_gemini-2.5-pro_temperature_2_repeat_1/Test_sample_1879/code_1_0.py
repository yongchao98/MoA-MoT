def solve_grid_paths():
    """
    Calculates the number of unique paths from (0,0) to (4,8) with specific
    move constraints using dynamic programming.
    """
    R = 4  # Target number of right moves
    U = 8  # Target number of up moves

    # dp[r][u] will store a list of 6 values:
    # [paths ending in R, RR, RRR, U, UU, UUU]
    dp = [[[0] * 6 for _ in range(U + 1)] for _ in range(R + 1)]

    # Base Cases:
    # Path to (1,0) is 'R' (ends in one R)
    if R > 0:
        dp[1][0][0] = 1
    # Path to (0,1) is 'U' (ends in one U)
    if U > 0:
        dp[0][1][3] = 1

    # Fill the DP table using the recurrence relations
    for r in range(R + 1):
        for u in range(U + 1):
            if r == 0 and u == 0:
                continue

            # Skip base cases that are already set
            if (r == 1 and u == 0) or (r == 0 and u == 1):
                continue

            # Calculate paths to (r, u)
            
            # Case 1: Path ends in one or more 'R's.
            # This requires coming from (r-1, u).
            if r > 0:
                # To end in a single R (...UR), the previous path at (r-1, u) must have ended in U.
                paths_from_u = dp[r-1][u][3] + dp[r-1][u][4] + dp[r-1][u][5]
                dp[r][u][0] = paths_from_u
                
                # To end in RR, the previous path at (r-1, u) must have ended in R.
                dp[r][u][1] = dp[r-1][u][0]
                
                # To end in RRR, the previous path at (r-1, u) must have ended in RR.
                dp[r][u][2] = dp[r-1][u][1]

            # Case 2: Path ends in one or more 'U's.
            # This requires coming from (r, u-1).
            if u > 0:
                # To end in a single U (...RU), the previous path at (r, u-1) must have ended in R.
                paths_from_r = dp[r][u-1][0] + dp[r][u-1][1] + dp[r][u-1][2]
                dp[r][u][3] = paths_from_r

                # To end in UU, the previous path at (r, u-1) must have ended in U.
                dp[r][u][4] = dp[r][u-1][3]
                
                # To end in UUU, the previous path at (r, u-1) must have ended in UU.
                dp[r][u][5] = dp[r][u-1][4]

    # The total number of paths is the sum of all possible valid endings at (R, U)
    final_counts = dp[R][U]
    total_paths = sum(final_counts)

    # Print the breakdown and the final result
    print("Number of paths ending with:")
    print(f"  ...R (1 R): {final_counts[0]}")
    print(f" ...RR (2 R's): {final_counts[1]}")
    print(f"...RRR (3 R's): {final_counts[2]}")
    print(f"  ...U (1 U): {final_counts[3]}")
    print(f" ...UU (2 U's): {final_counts[4]}")
    print(f"...UUU (3 U's): {final_counts[5]}")
    print("\nFinal Equation:")
    equation = f"{final_counts[0]} + {final_counts[1]} + {final_counts[2]} + {final_counts[3]} + {final_counts[4]} + {final_counts[5]} = {total_paths}"
    print(equation)
    print(f"\nTotal unique ways: {total_paths}")

solve_grid_paths()
<<<207>>>