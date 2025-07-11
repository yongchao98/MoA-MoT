def solve_coin_puzzle():
    """
    This function solves the coin puzzle and explains the strategy step-by-step.
    It calculates the maximum number of real coins that can be guaranteed to be identified.
    """
    total_coins = 1000
    fake_coins = 4

    # The optimal strategy uses specific group sizes.
    # We divide the 1000 coins into three groups: A, B, and C.
    size_A = 250
    size_B = 250
    size_C = total_coins - size_A - size_B

    # --- Weighing 1 Strategy ---
    # Weigh Group A vs Group B.

    # --- Analysis of Weighing 1 ---
    # Case 1: The scales balance (A = B).
    # All fakes must be in Group C. Groups A and B are entirely real.
    identified_if_balanced = size_A + size_B

    # Case 2: The scales are unbalanced (e.g., A < B).
    # This requires a second weighing. We split Group C into C1 and C2.
    size_C1 = 250
    size_C2 = 250

    # --- Weighing 2 Strategy (if A < B) ---
    # Weigh Group B vs Group C1.

    # --- Analysis of Weighing 2 ---
    # We find the guaranteed number of real coins for each outcome.
    # Let f_X be the number of fakes in group X. We know f_A > f_B.
    # Total fakes: f_A + f_B + f_C1 + f_C2 = 4.

    # Outcome 2a: B > C1 (B is heavier). This implies f_C1 > f_B.
    # A logical proof shows f_B must be 0. So, Group B is real.
    identified_if_B_gt_C1 = size_B

    # Outcome 2b: B < C1 (B is lighter). This implies f_B > f_C1.
    # A logical proof shows f_C1 must be 0. So, Group C1 is real.
    identified_if_B_lt_C1 = size_C1

    # Outcome 2c: B = C1 (they balance). This implies f_B = f_C1.
    # The worst-case distribution of fakes is (f_A, f_B, f_C1, f_C2) = (2, 1, 1, 0).
    # In this scenario, only Group C2 is guaranteed to be real.
    identified_if_B_eq_C1 = size_C2

    # The guaranteed number in the unbalanced case is the minimum of the three outcomes.
    guaranteed_in_unbalanced_case = min(identified_if_B_gt_C1, identified_if_B_lt_C1, identified_if_B_eq_C1)

    # The overall guarantee is the minimum of the two branches from Weighing 1.
    max_guaranteed_coins = min(identified_if_balanced, guaranteed_in_unbalanced_case)

    print("--- Coin Identification Strategy ---")
    print(f"We start with {total_coins} coins, including {fake_coins} lighter fake coins.")
    print("\nStep 1: Divide coins into three groups.")
    print(f"Group A: {size_A} coins")
    print(f"Group B: {size_B} coins")
    print(f"Group C: {size_C} coins")
    
    print("\n--- Weighing 1: Group A vs Group B ---")
    
    print("\nOutcome 1: The scales balance (A = B).")
    print("This implies all fake coins are in Group C.")
    print("Conclusion: All coins in Group A and Group B are real.")
    print(f"Number of identified real coins = {size_A} + {size_B} = {identified_if_balanced}")
    
    print("\nOutcome 2: The scales are unbalanced (e.g., A < B).")
    print("This means Group A has more fake coins than Group B. We proceed to Weighing 2.")
    print(f"We divide Group C ({size_C} coins) into C1 ({size_C1} coins) and C2 ({size_C2} coins).")
    
    print("\n--- Weighing 2 (if A < B): Group B vs Group C1 ---")
    
    print("We analyze all possible outcomes of this second weighing:")
    print(f"  - If B > C1 (B is heavier), it proves Group B must be real. We identify {identified_if_B_gt_C1} coins.")
    print(f"  - If B < C1 (B is lighter), it proves Group C1 must be real. We identify {identified_if_B_lt_C1} coins.")
    print(f"  - If B = C1 (they balance), a worst-case analysis shows that Group C2 must be real. We identify {identified_if_B_eq_C1} coins.")
    print(f"\nIn any outcome of Weighing 2, we are guaranteed to identify at least {guaranteed_in_unbalanced_case} real coins.")
    
    print("\n--- Final Conclusion ---")
    print("The number of real coins we can GUARANTEE to identify is the minimum of the results from the two main branches of the strategy.")
    print("Final guaranteed number = min(coins from balanced case, coins from unbalanced case)")
    print(f"Final guaranteed number = min({identified_if_balanced}, {guaranteed_in_unbalanced_case}) = {max_guaranteed_coins}")


solve_coin_puzzle()
<<<250>>>