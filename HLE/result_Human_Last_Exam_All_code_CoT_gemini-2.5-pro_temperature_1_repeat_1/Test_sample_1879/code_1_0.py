def solve_grid_walk():
    """
    Calculates the number of unique ways to walk from (0,0) to (4,8)
    with the constraint of no more than 3 consecutive steps in the same direction.
    """
    R_max = 4
    U_max = 8
    max_consecutive = 3

    # dp[r][u] will store {'R': [R1, R2, R3], 'U': [U1, U2, U3]}
    # Rx/Ux is the number of ways to get to (r,u) ending in x consecutive R/U moves.
    dp = [[{'R': [0] * max_consecutive, 'U': [0] * max_consecutive} for _ in range(U_max + 1)] for _ in range(R_max + 1)]

    # Base case: At (0,0), we assume a starting point that allows a first move
    # in either direction. A simpler way is to seed the first steps.
    if R_max > 0:
        dp[1][0]['R'][0] = 1  # Path "R" to (1,0)
    if U_max > 0:
        dp[0][1]['U'][0] = 1  # Path "U" to (0,1)

    # Iterate through the grid diagonally to ensure dependencies are met
    for s in range(2, R_max + U_max + 1):
        for r in range(s + 1):
            u = s - r
            if r > R_max or u > U_max:
                continue

            # Calculate ways to reach (r, u) ending with a Right move
            if r > 0:
                # To end with a single R, must come from a path ending in U
                total_from_U_at_prev = sum(dp[r - 1][u]['U'])
                dp[r][u]['R'][0] = total_from_U_at_prev
                # To end with RR or RRR, must come from a path ending in R or RR
                for k in range(1, max_consecutive):
                    dp[r][u]['R'][k] = dp[r - 1][u]['R'][k - 1]

            # Calculate ways to reach (r, u) ending with an Up move
            if u > 0:
                # To end with a single U, must come from a path ending in R
                total_from_R_at_prev = sum(dp[r][u - 1]['R'])
                dp[r][u]['U'][0] = total_from_R_at_prev
                # To end with UU or UUU, must come from a path ending in U or UU
                for k in range(1, max_consecutive):
                    dp[r][u]['U'][k] = dp[r][u - 1]['U'][k - 1]

    # The final answer is the sum of all ways to reach the destination (4,8)
    final_ways_R = dp[R_max][U_max]['R']
    final_ways_U = dp[R_max][U_max]['U']
    
    total_ways = sum(final_ways_R) + sum(final_ways_U)

    # Print the detailed breakdown of the final calculation
    r1, r2, r3 = final_ways_R
    u1, u2, u3 = final_ways_U
    
    print(f"Calculation Breakdown for paths to ({R_max},{U_max}):")
    print(f"Ways ending in 1 R: {r1}")
    print(f"Ways ending in 2 R's: {r2}")
    print(f"Ways ending in 3 R's: {r3}")
    print(f"Ways ending in 1 U: {u1}")
    print(f"Ways ending in 2 U's: {u2}")
    print(f"Ways ending in 3 U's: {u3}")
    print("\nFinal Equation:")
    print(f"{r1} + {r2} + {r3} + {u1} + {u2} + {u3} = {total_ways}")
    
    return total_ways

# Run the solver and capture the final answer
final_answer = solve_grid_walk()
print(f"\nTotal unique ways: {final_answer}")
print(f"<<<{final_answer}>>>")
