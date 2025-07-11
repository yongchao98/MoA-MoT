def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) with a maximum
    of 3 consecutive moves in the same direction.
    """
    WIDTH = 4
    HEIGHT = 8
    # "Cannot move four or more", so the maximum allowed consecutive moves is 3.
    MAX_CONSECUTIVE = 3

    # dp[x][y][dir][k]
    # x: x-coordinate (0 to WIDTH)
    # y: y-coordinate (0 to HEIGHT)
    # dir: 0 for Right, 1 for Up
    # k: number of consecutive steps (1 to MAX_CONSECUTIVE)
    dp = [[[[0] * (MAX_CONSECUTIVE + 1) for _ in range(2)] for _ in range(HEIGHT + 1)] for _ in range(WIDTH + 1)]

    # Base cases: paths along the axes
    # Paths along x-axis can only be sequences of 'R'
    for x in range(1, WIDTH + 1):
        if x <= MAX_CONSECUTIVE:
            dp[x][0][0][x] = 1
        else:
            break  # Any longer path on this axis is invalid

    # Paths along y-axis can only be sequences of 'U'
    for y in range(1, HEIGHT + 1):
        if y <= MAX_CONSECUTIVE:
            dp[0][y][1][y] = 1
        else:
            break  # Any longer path on this axis is invalid

    # Fill the rest of the DP table
    for x in range(1, WIDTH + 1):
        for y in range(1, HEIGHT + 1):
            # Calculate ways to reach (x, y) ending with Right moves (from (x-1, y))
            # Case 1: End with a single 'R'. This must follow a sequence of 'U's.
            paths_ending_up_at_prev = sum(dp[x-1][y][1][k] for k in range(1, MAX_CONSECUTIVE + 1))
            dp[x][y][0][1] = paths_ending_up_at_prev
            
            # Case 2: Extend a sequence of 'R's.
            for k in range(2, MAX_CONSECUTIVE + 1):
                dp[x][y][0][k] = dp[x-1][y][0][k-1]

            # Calculate ways to reach (x, y) ending with Up moves (from (x, y-1))
            # Case 1: End with a single 'U'. This must follow a sequence of 'R's.
            paths_ending_right_at_prev = sum(dp[x][y-1][0][k] for k in range(1, MAX_CONSECUTIVE + 1))
            dp[x][y][1][1] = paths_ending_right_at_prev

            # Case 2: Extend a sequence of 'U's.
            for k in range(2, MAX_CONSECUTIVE + 1):
                dp[x][y][1][k] = dp[x][y-1][1][k-1]

    # The final answer is the sum of all ways to get to (WIDTH, HEIGHT)
    final_counts = {}
    final_counts["ending in R"] = dp[WIDTH][HEIGHT][0][1]
    final_counts["ending in RR"] = dp[WIDTH][HEIGHT][0][2]
    final_counts["ending in RRR"] = dp[WIDTH][HEIGHT][0][3]
    final_counts["ending in U"] = dp[WIDTH][HEIGHT][1][1]
    final_counts["ending in UU"] = dp[WIDTH][HEIGHT][1][2]
    final_counts["ending in UUU"] = dp[WIDTH][HEIGHT][1][3]
    
    total_ways = 0
    sum_components = []
    
    print("The total number of unique ways is the sum of paths ending in different sequences:")
    for description, count in final_counts.items():
        if count > 0:
            print(f" - Paths {description}: {count}")
            total_ways += count
            sum_components.append(str(count))

    print("\nFinal calculation:")
    print(f"{' + '.join(sum_components)} = {total_ways}")
    
    return total_ways

if __name__ == '__main__':
    solve_grid_path()