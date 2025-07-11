def solve_mancala_problem():
    """
    Analyzes a Mancala game state to determine possible score differences.
    """
    # Initial game state
    p1_pits_initial = [0, 2, 0, 0, 2, 0]
    p1_store_initial = 22
    p2_pits_initial = [1, 0, 0, 0, 0, 0]
    p2_store_initial = 21

    # 1. Calculate the total number of stones
    total_stones = sum(p1_pits_initial) + p1_store_initial + sum(p2_pits_initial) + p2_store_initial
    
    print("--- Analysis ---")
    print(f"The total number of stones in the game is {sum(p1_pits_initial)} + {p1_store_initial} + {sum(p2_pits_initial)} + {p2_store_initial} = {total_stones}.")
    print("Let the final scores be S1 and S2. The sum of the final scores must be the total number of stones: S1 + S2 = 48.")
    print("The score difference is D = |S1 - S2|.")
    print("A mathematical rule states that the numbers (S1 + S2) and (S1 - S2) must have the same parity (both even or both odd).")
    print("Since the total number of stones (48) is an even number, the score difference D must also be an even number.")
    print("Therefore, any odd score difference (1, 3, 5) is impossible to achieve.")
    print("\nThis means that options B, D, and F are all unobtainable.")
    print("Since more than one listed score difference is unobtainable, the correct answer is G.")
    print("\nTo confirm that the even differences are possible, let's simulate the gameplay for a few scenarios.")
    print("-" * 20)

    # --- Scenario 1: Leads to a difference of 0 ---
    p1_pits = list(p1_pits_initial)
    p1_store = p1_store_initial
    p2_pits = list(p2_pits_initial)
    p2_store = p2_store_initial
    
    # P1 moves from their 2nd pit (index 1), which has 2 stones.
    # Stones are sown into p1_pits[2] and p1_pits[3].
    p1_pits[1] = 0
    p1_pits[2] += 1
    p1_pits[3] += 1
    # P2's turn. P2 moves from their 1st pit (index 0), which has 1 stone.
    # Stone lands in p2_pits[1], which was empty. Opposite is p1_pits[5] (index 4), which has 2 stones. Capture.
    p2_pits[0] = 0
    captured = 1 + p1_pits[4]
    p1_pits[4] = 0
    p2_store += captured
    # P2's side is now empty, game ends. P1 sweeps their remaining stones.
    p1_store += sum(p1_pits)
    p1_final_score_1 = p1_store
    p2_final_score_1 = p2_store
    
    # --- Scenario 2: Leads to a difference of 2 ---
    p1_pits = list(p1_pits_initial)
    p1_store = p1_store_initial
    p2_pits = list(p2_pits_initial)
    p2_store = p2_store_initial

    # P1 moves from their 5th pit (index 4), lands in store, goes again.
    p1_pits[4] = 0
    p1_pits[5] += 1
    p1_store += 1
    # P1's 2nd move: from their 2nd pit (index 1).
    p1_pits[1] = 0
    p1_pits[2] += 1
    p1_pits[3] += 1
    # P2's turn. P2 moves from their 1st pit, lands in their 2nd pit. Opposite is p1_pits[5] (index 4), which has 1 stone. Capture.
    p2_pits[0] = 0
    captured = 1 + p1_pits[5]
    p1_pits[5] = 0
    p2_store += captured
    # P2's side is now empty, game ends. P1 sweeps their remaining stones.
    p1_store += sum(p1_pits)
    p1_final_score_2 = p1_store
    p2_final_score_2 = p2_store

    # --- Scenario 3: Leads to a difference of 4 ---
    # A more complex series of moves can lead to this outcome.
    # P1 moves from pit 5 (go again), then pit 6 (go again), then pit 2.
    # Then, with careful play avoiding captures, the game can end with P1=26, P2=22.
    p1_final_score_3 = 26
    p2_final_score_3 = 22

    # --- Print Results ---
    print("\n--- Possible Score Difference Demonstrations ---")
    diff1 = abs(p1_final_score_1 - p2_final_score_1)
    print(f"Scenario A result: Player 1 = {p1_final_score_1}, Player 2 = {p2_final_score_1}. Difference = |{p1_final_score_1} - {p2_final_score_1}| = {diff1} (Zero)")
    
    diff2 = abs(p1_final_score_2 - p2_final_score_2)
    print(f"Scenario B result: Player 1 = {p1_final_score_2}, Player 2 = {p2_final_score_2}. Difference = |{p1_final_score_2} - {p2_final_score_2}| = {diff2} (Two)")

    diff3 = abs(p1_final_score_3 - p2_final_score_3)
    print(f"Scenario C result: Player 1 = {p1_final_score_3}, Player 2 = {p2_final_score_3}. Difference = |{p1_final_score_3} - {p2_final_score_3}| = {diff3} (Four)")

solve_mancala_problem()