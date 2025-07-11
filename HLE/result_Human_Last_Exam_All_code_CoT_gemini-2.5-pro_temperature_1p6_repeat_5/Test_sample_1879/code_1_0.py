def solve_grid_path():
    """
    Calculates the number of unique paths from (0,0) to (4,8) on a 2D grid
    with the constraint of no more than 3 consecutive moves in the same direction.
    """
    x_max, y_max = 4, 8
    k_max = 3  # Max consecutive moves allowed

    # dp[x][y][dir][k]
    # x: x-coordinate (0 to x_max)
    # y: y-coordinate (0 to y_max)
    # dir: 0 for Right, 1 for Up
    # k: consecutive moves in dir (1 to k_max)
    # Initialize DP table with zeros
    dp = [[[[0] * (k_max + 1) for _ in range(2)] for _ in range(y_max + 1)] for _ in range(x_max + 1)]

    # Iterate through the grid, diagonal by diagonal, ensuring dependencies are met
    for s in range(1, x_max + y_max + 1):
        for x in range(x_max + 1):
            y = s - x
            if not (0 <= y <= y_max):
                continue

            # Base cases: First move from origin (0,0)
            if x == 1 and y == 0:
                dp[1][0][0][1] = 1
                continue
            if x == 0 and y == 1:
                dp[0][1][1][1] = 1
                continue

            # Calculate ways ending in RIGHT moves, arriving at (x,y) from (x-1,y)
            if x > 0:
                # Arriving with 1 Right move (sequence ends like ...UR)
                ways_from_up = 0
                for k in range(1, k_max + 1):
                    ways_from_up += dp[x-1][y][1][k]
                dp[x][y][0][1] = ways_from_up
                
                # Arriving with 2 or more Right moves (sequence ends like ...RR, ...RRR)
                for k in range(2, k_max + 1):
                    dp[x][y][0][k] = dp[x-1][y][0][k-1]

            # Calculate ways ending in UP moves, arriving at (x,y) from (x,y-1)
            if y > 0:
                # Arriving with 1 Up move (sequence ends like ...RU)
                ways_from_right = 0
                for k in range(1, k_max + 1):
                    ways_from_right += dp[x][y-1][0][k]
                dp[x][y][1][1] = ways_from_right

                # Arriving with 2 or more Up moves (sequence ends like ...UU, ...UUU)
                for k in range(2, k_max + 1):
                    dp[x][y][1][k] = dp[x][y-1][1][k-1]

    # The final answer is the sum of all ways to arrive at the destination (x_max, y_max)
    w = {
        'r1': dp[x_max][y_max][0][1],
        'r2': dp[x_max][y_max][0][2],
        'r3': dp[x_max][y_max][0][3],
        'u1': dp[x_max][y_max][1][1],
        'u2': dp[x_max][y_max][1][2],
        'u3': dp[x_max][y_max][1][3],
    }
    
    total_ways = sum(w.values())

    print("The total number of unique ways is the sum of ways to arrive at (4,8):")
    print(f" - Ways ending in 1 Right move (...UR): {w['r1']}")
    print(f" - Ways ending in 2 Right moves (...RR): {w['r2']}")
    print(f" - Ways ending in 3 Right moves (...RRR): {w['r3']}")
    print(f" - Ways ending in 1 Up move (...RU): {w['u1']}")
    print(f" - Ways ending in 2 Up moves (...UU): {w['u2']}")
    print(f" - Ways ending in 3 Up moves (...UUU): {w['u3']}")
    print("\nFinal Equation:")
    print(f"{w['r1']} + {w['r2']} + {w['r3']} + {w['u1']} + {w['u2']} + {w['u3']} = {total_ways}")
    print(f"\n<<< {total_ways} >>>")

solve_grid_path()