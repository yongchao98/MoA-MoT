import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at the conference table
    based on a complex set of rules.
    """
    factorial = math.factorial

    # Step 1: Define group sizes and compositions
    n_scientists = 12
    n_mathematicians = 4
    n_ethicists = 2
    n_classicists = 6  # 5 picky classicists + Cassie
    n_picky_classicists = 5

    w_scientists = n_scientists // 2  # 6
    m_scientists = n_scientists // 2  # 6
    w_mathematicians = n_mathematicians // 2  # 2
    m_mathematicians = n_mathematicians // 2  # 2

    # Rowers are all male
    rower_s = 2
    rower_m = 1

    # Calculate remaining "other" members in the Scientist and Mathematician groups
    other_s = n_scientists - rower_s  # 10
    other_w_s = w_scientists  # 6
    other_m_s = m_scientists - rower_s  # 4

    other_m = n_mathematicians - rower_m  # 3
    other_w_m = w_mathematicians  # 2
    other_m_m = m_mathematicians - rower_m  # 1

    # Step 2: Calculate permutations for internal block components
    # The rower block {S,S,M} must be together. The two scientists {S,S} form a sub-block.
    # Ways to arrange the two sub-units {(SS), M} is 2!. Ways to arrange scientists within their sub-block is 2!
    ways_R = factorial(2) * factorial(2)

    # Calculate permutations for the "other" groups, considering the gender of the person at the end of the line
    # "Other Scientists" (10 people: 6W, 4M)
    ways_S_other_total = factorial(other_s)
    ways_S_other_female_end = other_w_s * factorial(other_s - 1)
    
    # "Other Mathematicians" (3 people: 2W, 1M)
    ways_M_other_total = factorial(other_m)
    ways_M_other_female_end = other_w_m * factorial(other_m - 1)

    # Step 3: Calculate permutations for external arrangements
    # "Gap" people are the 5 picky classicists + 1 leftover buffer person
    n_gap_people = n_picky_classicists + 1
    ways_gap = factorial(n_gap_people)

    # Step 4: Sum the arrangements for each valid case
    
    # Case 1: Neighbors of the SM_super_block are the two Ethicists.
    # External ways: 2! for arranging ethicists * 6! for the gap.
    ext_ways_E_E = factorial(n_ethicists) * ways_gap
    # Internal ways: No gender constraint on the ends of the block.
    int_ways_E_E = ways_R * ways_S_other_total * ways_M_other_total
    total_case1 = ext_ways_E_E * int_ways_E_E

    # Case 2: Neighbors are Cassie and an Ethicist, with Cassie at the Scientist-end.
    # S-end must be female.
    # External ways: Choose 1 of 2 ethicists * 6! for the gap.
    ext_ways_K_A_E = n_ethicists * ways_gap
    # Internal ways: S-end must be female.
    int_ways_K_A_E = ways_R * ways_S_other_female_end * ways_M_other_total
    total_case2 = ext_ways_K_A_E * int_ways_K_A_E

    # Case 3: Neighbors are an Ethicist and Cassie, with Cassie at the Mathematician-end.
    # M-end must be female.
    # External ways: Choose 1 of 2 ethicists * 6! for the gap.
    ext_ways_E_A_K = n_ethicists * ways_gap
    # Internal ways: M-end must be female.
    int_ways_E_A_K = ways_R * ways_S_other_total * ways_M_other_female_end
    total_case3 = ext_ways_E_A_K * int_ways_E_A_K
    
    # Final Result
    total_ways = total_case1 + total_case2 + total_case3

    # Print the breakdown of the calculation
    print("The total number of arrangements is the sum of three cases:")
    print("-" * 50)
    
    # Print Case 1
    f_ext_E_E = f"({factorial(n_ethicists)} * {ways_gap})"
    f_int_E_E = f"({ways_R} * {ways_S_other_total} * {ways_M_other_total})"
    print(f"1. Ethicist Neighbors: {f_ext_E_E} * {f_int_E_E} = {total_case1}")
    
    # Print Case 2
    f_ext_K_A_E = f"({n_ethicists} * {ways_gap})"
    f_int_K_A_E = f"({ways_R} * ({other_w_s} * {factorial(other_s - 1)}) * {ways_M_other_total})"
    print(f"2. Cassie at Scientist End: {f_ext_K_A_E} * {f_int_K_A_E} = {total_case2}")

    # Print Case 3
    f_ext_E_A_K = f"({n_ethicists} * {ways_gap})"
    f_int_E_A_K = f"({ways_R} * {ways_S_other_total} * ({other_w_m} * {factorial(other_m - 1)}))"
    print(f"3. Cassie at Mathematician End: {f_ext_E_A_K} * {f_int_E_A_K} = {total_case3}")

    print("-" * 50)
    print(f"Total ways = {total_case1} + {total_case2} + {total_case3}")
    print(f"Total ways = {total_ways}")

solve_seating_arrangement()
<<<284265676800>>>