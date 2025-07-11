import collections

def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    R_MAX = 4
    U_MAX = 8
    CONSECUTIVE_MAX = 3

    # dp[r][u][dir][len]
    # dir: 0 for Right, 1 for Up
    # len: 1, 2, or 3 for the length of the consecutive moves.
    # We use a table of size CONSECUTIVE_MAX + 1 for 1-based indexing of len.
    dp = [[[[0] * (CONSECUTIVE_MAX + 1) for _ in range(2)] for _ in range(U_MAX + 1)] for _ in range(R_MAX + 1)]

    # Base cases for the first move from (0,0)
    if R_MAX > 0:
        # One way to get to (1,0) is by one 'Right' move.
        dp[1][0][0][1] = 1
    if U_MAX > 0:
        # One way to get to (0,1) is by one 'Up' move.
        dp[0][1][1][1] = 1

    # Fill the DP table by iterating through each cell of the grid
    for r in range(R_MAX + 1):
        for u in range(U_MAX + 1):
            if r == 0 and u == 0:
                continue

            # Calculate ways to reach (r, u) ending with a Right move.
            # This must come from cell (r-1, u).
            if r > 0:
                # Case 1: The last move to (r-1, u) was Up. We start a new sequence of Right moves (length 1).
                # Sum all ways to get to (r-1, u) that ended in an Up move.
                ways_from_up = sum(dp[r - 1][u][1])
                dp[r][u][0][1] = ways_from_up
                
                # Case 2: The last move to (r-1, u) was also Right. We extend the sequence.
                for length in range(2, CONSECUTIVE_MAX + 1):
                    dp[r][u][0][length] = dp[r-1][u][0][length - 1]

            # Calculate ways to reach (r, u) ending with an Up move.
            # This must come from cell (r, u-1).
            if u > 0:
                # Case 1: The last move to (r, u-1) was Right. We start a new sequence of Up moves (length 1).
                # Sum all ways to get to (r, u-1) that ended in a Right move.
                ways_from_right = sum(dp[r][u - 1][0])
                dp[r][u][1][1] = ways_from_right

                # Case 2: The last move to (r, u-1) was also Up. We extend the sequence.
                for length in range(2, CONSECUTIVE_MAX + 1):
                    dp[r][u][1][length] = dp[r][u-1][1][length - 1]

    # The final answer is the sum of all ways to reach the target cell (R_MAX, U_MAX).
    final_ways_ending_in_R = dp[R_MAX][U_MAX][0]
    final_ways_ending_in_U = dp[R_MAX][U_MAX][1]
    
    r1, r2, r3 = final_ways_ending_in_R[1], final_ways_ending_in_R[2], final_ways_ending_in_R[3]
    u1, u2, u3 = final_ways_ending_in_U[1], final_ways_ending_in_U[2], final_ways_ending_in_U[3]

    total_ways = r1 + r2 + r3 + u1 + u2 + u3

    print(f"Paths to ({R_MAX},{U_MAX}) ending with 1 Right move: {r1}")
    print(f"Paths to ({R_MAX},{U_MAX}) ending with 2 Right moves: {r2}")
    print(f"Paths to ({R_MAX},{U_MAX}) ending with 3 Right moves: {r3}")
    print(f"Paths to ({R_MAX},{U_MAX}) ending with 1 Up move: {u1}")
    print(f"Paths to ({R_MAX},{U_MAX}) ending with 2 Up moves: {u2}")
    print(f"Paths to ({R_MAX},{U_MAX}) ending with 3 Up moves: {u3}")
    print(f"Total unique ways = {r1} + {r2} + {r3} + {u1} + {u2} + {u3} = {total_ways}")
    
    return total_ways

final_answer = solve_grid_path()
print(f"<<<{final_answer}>>>")