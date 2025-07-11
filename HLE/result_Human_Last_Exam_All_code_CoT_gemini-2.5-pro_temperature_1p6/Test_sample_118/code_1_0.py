def solve_coin_problem():
    """
    Calculates the maximum number of real coins guaranteed to be identified.
    This function implements the logic described above.
    """
    
    total_coins = 1000
    fake_coins = 4

    # Weighing 1: 1 coin (A) vs 1 coin (B). Remainder (C) is 998 coins.
    
    # Case 1: First weighing is unequal (e.g., A < B)
    # This means A is fake, B is real. We've identified 1 real coin.
    # The remaining 3 fakes are in the 998 coins of group C.
    real_coins_from_w1_unequal = 1
    coins_in_c = 998
    fakes_in_c_unequal = 3
    
    # Weighing 2 for this case: Split C into 333, 333, 332.
    # The worst case for identifying real coins is when the suspect group is largest.
    # If we weigh 333 vs 333, the fakes could be in one group of 333.
    c_group1_size = 333
    c_group2_size = 333
    c_group3_size = 332
    worst_suspect_group_unequal = max(c_group1_size, c_group2_size, c_group3_size)
    
    # The minimum number of real coins identified in group C
    identified_in_c_unequal = coins_in_c - worst_suspect_group_unequal
    
    # Total guaranteed real coins in this path
    total_guaranteed_unequal = real_coins_from_w1_unequal + identified_in_c_unequal
    
    print("--- Analysis for Unequal First Weighing (e.g. A < B) ---")
    print(f"Step 1: Weigh 1 coin vs 1 coin. The heavier coin is real.")
    print(f"Real coins identified: {real_coins_from_w1_unequal}")
    print(f"The remaining {coins_in_c} coins contain {fakes_in_c_unequal} fakes.")
    print(f"Step 2: Weigh two groups of {c_group1_size} from the remainder against each other.")
    print(f"Worst case: all fakes are in one of these groups, making {worst_suspect_group_unequal} coins suspect.")
    print(f"Guaranteed real coins from remainder = {coins_in_c} - {worst_suspect_group_unequal} = {identified_in_c_unequal}")
    print(f"Total guaranteed real coins if W1 is unequal = {real_coins_from_w1_unequal} + {identified_in_c_unequal} = {total_guaranteed_unequal}")
    print("-" * 20)
    
    # Case 2: First weighing is equal (A = B)
    # This means A and B are both real. We've identified 2 real coins.
    # The 4 fakes are in the 998 coins of group C.
    real_coins_from_w1_equal = 2
    fakes_in_c_equal = 4
    
    # Weighing 2 for this case is the same: C1(333) vs C2(333).
    # The worst case again is that the suspect group is the largest possible.
    worst_suspect_group_equal = max(c_group1_size, c_group2_size)
    identified_in_c_equal = coins_in_c - worst_suspect_group_equal
    
    # Total guaranteed real coins in this path
    total_guaranteed_equal = real_coins_from_w1_equal + identified_in_c_equal
    
    print("--- Analysis for Equal First Weighing (A = B) ---")
    print(f"Step 1: Weigh 1 coin vs 1 coin. If they balance, both are real.")
    print(f"Real coins identified: {real_coins_from_w1_equal}")
    print(f"The remaining {coins_in_c} coins contain {fakes_in_c_equal} fakes.")
    print(f"Step 2: Weigh two groups of {c_group1_size} from the remainder against each other.")
    print(f"Worst case: all fakes are in one of these groups, making {worst_suspect_group_equal} coins suspect.")
    print(f"Guaranteed real coins from remainder = {coins_in_c} - {worst_suspect_group_equal} = {identified_in_c_equal}")
    print(f"Total guaranteed real coins if W1 is equal = {real_coins_from_w1_equal} + {identified_in_c_equal} = {total_guaranteed_equal}")
    print("-" * 20)

    # The number of coins we can *guarantee* to identify is the minimum of the outcomes.
    max_guaranteed_reals = min(total_guaranteed_unequal, total_guaranteed_equal)
    
    print("The maximum number of real coins we can GUARANTEE to identify is the minimum of the results from all possible outcomes.")
    print(f"Final Answer = min({total_guaranteed_unequal}, {total_guaranteed_equal}) = {max_guaranteed_reals}")
    
    return max_guaranteed_reals

solve_coin_problem()