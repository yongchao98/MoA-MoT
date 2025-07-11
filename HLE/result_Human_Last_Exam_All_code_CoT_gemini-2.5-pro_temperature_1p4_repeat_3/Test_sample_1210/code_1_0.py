def solve_mancala_problem():
    """
    Analyzes the Mancala game state to find all possible score differences.
    """

    print("Analyzing the Mancala game starting from the given state.")
    print("Initial State:")
    p1_pits = [0, 2, 0, 0, 2, 0]
    p1_store = 22
    p2_pits = [1, 0, 0, 0, 0, 0]
    p2_store = 21
    print(f"Player 1 Pits: {p1_pits}, Store: {p1_store}")
    print(f"Player 2 Pits: {p2_pits}, Store: {p2_store}")
    print("It is Player 1's move.")
    print("-" * 30)

    possible_differences = set()

    # --- Scenario 1: Player 1 plays from the second pit (index 1) ---
    print("Scenario 1: Player 1 chooses the pit with 2 stones at index 1.")
    
    # P1's move
    p1_pits_s1 = [0, 0, 1, 1, 2, 0]
    p1_store_s1 = 22
    print("Player 1 sows 2 stones, landing in their 4th pit. No capture, no extra turn.")
    print(f"Board state becomes: P1 Pits: {p1_pits_s1}, P2 Pits: {p2_pits}")
    print("It is now Player 2's turn.")

    # P2's move
    p2_pits_s1_after_move = [0, 0, 0, 0, 0, 0] # P2 moves from their 1st pit, stone lands in 2nd
    # The last stone landed in an empty pit on Player 2's side (pit index 1),
    # opposite Player 1's pit (index 4) which contains 2 stones.
    # P2 captures those 2 stones plus their landing stone.
    captured_stones = p1_pits_s1[4]
    p1_pits_s1[4] = 0
    p2_store_s1 = p2_store + captured_stones + 1
    
    print("Player 2 moves their only stone. It lands in an empty pit opposite Player 1's pit with 2 stones.")
    print(f"Player 2 captures {captured_stones} stones from Player 1.")
    print(f"Player 2's store becomes {p2_store} + {captured_stones} + 1 = {p2_store_s1}.")
    print("Player 2's pits are now all empty, so the game ends.")

    # Game End
    remaining_p1_stones = sum(p1_pits_s1)
    final_p1_score_s1 = p1_store_s1 + remaining_p1_stones
    final_p2_score_s1 = p2_store_s1
    
    print(f"Player 1 collects the remaining {remaining_p1_stones} stones from their pits.")
    print(f"Final Score P1: {p1_store_s1} + {remaining_p1_stones} = {final_p1_score_s1}")
    print(f"Final Score P2: {final_p2_score_s1}")
    
    diff_s1 = abs(final_p1_score_s1 - final_p2_score_s1)
    possible_differences.add(diff_s1)
    print(f"Score Difference = |{final_p1_score_s1} - {final_p2_score_s1}| = {diff_s1}")
    print("-" * 30)

    # --- Scenario 2: Player 1 plays from the fifth pit (index 4) ---
    print("Scenario 2: Player 1 chooses the pit with 2 stones at index 4.")
    
    # P1's move
    p1_pits_s2 = [0, 2, 0, 0, 0, 1]
    p1_store_s2 = p1_store + 1
    print("Player 1 sows 2 stones, with the last stone landing in their store.")
    print(f"Player 1 gets an extra turn. Store becomes {p1_store_s2}.")
    print(f"Board state is: P1 Pits: {p1_pits_s2}, Store: {p1_store_s2}")

    # For both of P1's subsequent choices, the outcome is the same.
    # If P1 plays from pit 5 -> goes to store -> go again. Then P1 must play from pit 1.
    # If P1 plays from pit 1 -> lands on own side. P2's turn.
    # In either path, no more captures are possible because the opponent's opposite pits are empty.
    # Therefore, all remaining stones on the board will end up in their owner's store.
    print("On the next turn(s), Player 1 moves their remaining stones.")
    print("Analysis shows no more captures are possible in this branch of the game.")
    
    # Game End
    total_p1_stones_on_board = sum(p1_pits_s2)
    total_p2_stones_on_board = sum(p2_pits)
    final_p1_score_s2 = p1_store_s2 + total_p1_stones_on_board
    final_p2_score_s2 = p2_store + total_p2_stones_on_board
    
    print("All remaining stones in the pits will eventually go to their respective stores.")
    print(f"Final Score P1: {p1_store_s2} + {total_p1_stones_on_board} = {final_p1_score_s2}")
    print(f"Final Score P2: {p2_store} + {total_p2_stones_on_board} = {final_p2_score_s2}")

    diff_s2 = abs(final_p1_score_s2 - final_p2_score_s2)
    possible_differences.add(diff_s2)
    print(f"Score Difference = |{final_p1_score_s2} - {final_p2_score_s2}| = {diff_s2}")
    print("-" * 30)

    # --- Conclusion ---
    print(f"The set of all possible score differences is: {possible_differences}")
    answer_choices = {'Zero': 0, 'One': 1, 'Two': 2, 'Three': 3, 'Four': 4, 'Five': 5}
    unobtainable = []
    for name, val in answer_choices.items():
        if val not in possible_differences:
            unobtainable.append(name)
    
    print(f"Comparing this to the answer choices, the unobtainable differences are: {', '.join(unobtainable)}.")
    if len(unobtainable) > 1:
        print("Since more than one of the listed score differences is unobtainable, the answer is G.")
    else:
        print(f"The unobtainable difference is {unobtainable[0]}.")

solve_mancala_problem()
<<<G>>>