def solve_coin_game():
    """
    Calculates the outcome of the coin game based on Player 1's optimal strategy.
    """
    # Initial game setup
    n1_coins = 136
    n2_coins = 87
    
    # Total turns per player
    total_coins = n1_coins + n2_coins
    p1_turns = (total_coins + 1) // 2
    p2_turns = total_coins // 2

    print("--- Game Analysis ---")
    print(f"Total coins: {total_coins}. P1 picks {p1_turns}, P2 picks {p2_turns}.")
    print("The score can be simplified:")
    print(f"Player 1 Score = {p1_turns} + (number of 2-euro coins P1 gets)")
    print(f"Player 2 Score = {p2_turns} + (number of 2-euro coins P2 gets)")
    print("The game is a race to collect the most 2-euro coins.")
    print("-" * 25)

    # --- Scenario 1: P1 takes a 2-euro coin first ---
    print("Scenario 1: Player 1's first move is to take a 2-euro coin.")
    n2_subgame_1 = n2_coins - 1
    # P1 gets 1 coin, P2 gets 0. Subgame starts.
    # In the even-item subgame, the items are split.
    p1_n2_sub_1 = n2_subgame_1 // 2
    p2_n2_sub_1 = n2_subgame_1 // 2
    
    p1_total_n2_1 = 1 + p1_n2_sub_1
    p2_total_n2_1 = p2_n2_sub_1
    
    p1_score_1 = p1_turns + p1_total_n2_1
    p2_score_1 = p2_turns + p2_total_n2_1

    print(f"  P1's total 2-euro coins: 1 (initial) + {p1_n2_sub_1} (subgame) = {p1_total_n2_1}")
    print(f"  P2's total 2-euro coins: 0 (initial) + {p2_n2_sub_1} (subgame) = {p2_total_n2_1}")
    print(f"  P1's final score: {p1_turns} + {p1_total_n2_1} = {p1_score_1}")
    print(f"  P2's final score: {p2_turns} + {p2_total_n2_1} = {p2_score_1}")
    print("  Outcome: Player 1 wins.")
    print("-" * 25)

    # --- Scenario 2: P1 must take a 1-euro coin first ---
    print("Scenario 2: Player 1's first move is to take a 1-euro coin.")
    n2_subgame_2 = n2_coins
    # In the odd-item subgame, the starter (P2) gets the extra coin.
    p1_n2_sub_2 = n2_subgame_2 // 2
    p2_n2_sub_2 = (n2_subgame_2 + 1) // 2

    p1_total_n2_2 = 0 + p1_n2_sub_2
    p2_total_n2_2 = p2_n2_sub_2

    p1_score_2 = p1_turns + p1_total_n2_2
    p2_score_2 = p2_turns + p2_total_n2_2
    
    print(f"  P1's total 2-euro coins: 0 (initial) + {p1_n2_sub_2} (subgame) = {p1_total_n2_2}")
    print(f"  P2's total 2-euro coins: 0 (initial) + {p2_n2_sub_2} (subgame) = {p2_total_n2_2}")
    print(f"  P1's final score: {p1_turns} + {p1_total_n2_2} = {p1_score_2}")
    print(f"  P2's final score: {p2_turns} + {p2_total_n2_2} = {p2_score_2}")
    print("  Outcome: Tie.")
    print("-" * 25)
    
    print("\nConclusion: Player 1 can always guarantee at least a tie, and will win if a 2-euro coin is initially at either end.")
    print("Therefore, you should choose to be the 1st player.")

solve_coin_game()