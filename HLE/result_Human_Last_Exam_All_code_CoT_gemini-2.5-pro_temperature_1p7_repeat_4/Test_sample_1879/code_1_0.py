def solve_grid_path():
    """
    Calculates the number of unique paths on a 2D grid from (0,0) to (4,8)
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    target_x = 4
    target_y = 8
    # The maximum number of consecutive moves allowed in the same direction.
    max_consecutive = 3

    # DP State: dp[x][y][dir][count]
    # x: x-coordinate (0 to 4)
    # y: y-coordinate (0 to 8)
    # dir: 0 for Right, 1 for Up
    # count: Number of consecutive moves (1 to 3)
    # We use size max_consecutive + 1 for 1-based indexing of count.
    dp = [[[[0] * (max_consecutive + 1) for _ in range(2)] for _ in range(target_y + 1)] for _ in range(target_x + 1)]

    # Iterate through all grid points to calculate the number of ways to ARRIVE at each one.
    for x in range(target_x + 1):
        for y in range(target_y + 1):
            if x == 0 and y == 0:
                continue

            # --- Calculate ways to arrive at (x,y) via a RIGHT move (from the left) ---
            if x > 0:
                prev_x, prev_y = x - 1, y

                # Case 1: Start a new sequence of Right moves.
                # The previous move to (prev_x, prev_y) must have been Up.
                # A path from (0,0) to (1,0) is a special base case.
                ways_from_left_ending_in_up = 0
                if prev_x == 0 and prev_y == 0:
                    ways_from_left_ending_in_up = 1
                else:
                    for k in range(1, max_consecutive + 1):
                        ways_from_left_ending_in_up += dp[prev_x][prev_y][1][k]
                dp[x][y][0][1] = ways_from_left_ending_in_up

                # Case 2: Continue a sequence of Right moves.
                for k in range(2, max_consecutive + 1):
                    dp[x][y][0][k] = dp[prev_x][prev_y][0][k - 1]

            # --- Calculate ways to arrive at (x,y) via an UP move (from below) ---
            if y > 0:
                prev_x, prev_y = x, y - 1

                # Case 1: Start a new sequence of Up moves.
                # The previous move to (prev_x, prev_y) must have been Right.
                # A path from (0,0) to (0,1) is a special base case.
                ways_from_down_ending_in_right = 0
                if prev_x == 0 and prev_y == 0:
                    ways_from_down_ending_in_right = 1
                else:
                    for k in range(1, max_consecutive + 1):
                        ways_from_down_ending_in_right += dp[prev_x][prev_y][0][k]
                dp[x][y][1][1] = ways_from_down_ending_in_right
                
                # Case 2: Continue a sequence of Up moves.
                for k in range(2, max_consecutive + 1):
                    dp[x][y][1][k] = dp[prev_x][prev_y][1][k - 1]
    
    # The final answer is the sum of all possible states at the target cell.
    sum_components = []
    total_ways = 0
    for direction in range(2):
        for count in range(1, max_consecutive + 1):
            ways = dp[target_x][target_y][direction][count]
            total_ways += ways
            sum_components.append(ways)

    # Print the detailed breakdown of the final sum
    equation_str = " + ".join(map(str, sum_components))
    print("The final answer is the sum of ways to arrive at (4,8) ending with:")
    print(f"- 1, 2, or 3 Right moves: {sum_components[0]}, {sum_components[1]}, {sum_components[2]}")
    print(f"- 1, 2, or 3 Up moves:    {sum_components[3]}, {sum_components[4]}, {sum_components[5]}")
    print("\nFinal calculation:")
    print(f"{equation_str} = {total_ways}")
    
    print("\nThe total number of unique ways is:")
    print(total_ways)
    return total_ways

# Execute the function to find the solution.
solve_grid_path()