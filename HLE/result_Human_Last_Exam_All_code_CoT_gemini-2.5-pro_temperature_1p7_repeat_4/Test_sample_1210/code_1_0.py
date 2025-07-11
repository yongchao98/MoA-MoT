import collections

def solve_mancala_puzzle():
    """
    Solves the Mancala puzzle by exploring all possible game states
    from the initial configuration to find all possible score differences.
    """

    # ------------------
    # 1. INITIAL SETUP
    # ------------------
    # The board is represented as a list of 14 elements.
    # P1 Pits: 0-5, P1 Store: 6
    # P2 Pits: 7-12, P2 Store: 13
    initial_board = [0, 2, 0, 0, 2, 0, 22, 1, 0, 0, 0, 0, 0, 21]
    p1_pits_indices = range(0, 6)
    p1_store_index = 6
    p2_pits_indices = range(7, 13)
    p2_store_index = 13
    
    # ------------------
    # 2. LOGICAL REASONING (PARITY ARGUMENT)
    # ------------------
    total_stones = sum(initial_board)
    print("Step 1: Logical Analysis\n")
    print(f"The initial state is P1 pits: {initial_board[0:6]}, P1 store: {initial_board[6]}; P2 pits: {initial_board[7:13]}, P2 store: {initial_board[13]}.")
    print(f"The total number of stones in the game is {sum(initial_board[0:7])} + {sum(initial_board[7:14])} = {total_stones}.\n")
    print("Let the final scores for Player 1 and Player 2 be S1 and S2.")
    print(f"At the end of the game, all {total_stones} stones will be in the stores, so S1 + S2 = {total_stones}.")
    print("The score difference is D = |S1 - S2|.")
    print("If we assume S1 is the winner's score, then S1 - S2 = D.")
    print("Adding the two equations (S1 + S2) + (S1 - S2) = 48 + D, we get 2*S1 = 48 + D.")
    print("This means S1 = (48 + D) / 2.\n")
    print("Since a player's score (S1) must be a whole number, (48 + D) must be an even number.")
    print("Because 48 is even, the score difference D must also be an EVEN number.")
    print("This mathematical property means that any odd score difference (1, 3, 5) is impossible to achieve.")
    print("This strongly suggests that more than one of the listed choices are unobtainable.\n")

    # ------------------
    # 3. GAME SIMULATION
    # ------------------
    print("Step 2: Game Simulation to Confirm\n")
    print("To verify and check the even-numbered choices, we will simulate all possible games from the start.")
    
    final_score_differences = set()
    queue = collections.deque([(tuple(initial_board), 1)])  # (board_tuple, player_to_move)
    visited = set([(tuple(initial_board), 1)])

    while queue:
        current_board_tuple, current_player = queue.popleft()
        current_board = list(current_board_tuple)

        # Check for game over
        p1_pits_empty = all(current_board[i] == 0 for i in p1_pits_indices)
        p2_pits_empty = all(current_board[i] == 0 for i in p2_pits_indices)

        if p1_pits_empty or p2_pits_empty:
            final_board = list(current_board)
            for i in p1_pits_indices:
                final_board[p1_store_index] += final_board[i]; final_board[i] = 0
            for i in p2_pits_indices:
                final_board[p2_store_index] += final_board[i]; final_board[i] = 0
            
            p1_score = final_board[p1_store_index]
            p2_score = final_board[p2_store_index]
            final_score_differences.add(abs(p1_score - p2_score))
            continue

        # Get and perform possible moves
        possible_move_pits = [i for i in (p1_pits_indices if current_player == 1 else p2_pits_indices) if current_board[i] > 0]
            
        for move_pit in possible_move_pits:
            board_after_move = list(current_board)
            stones = board_after_move[move_pit]
            board_after_move[move_pit] = 0

            # Distribute stones
            last_pit_index = move_pit
            for _ in range(stones):
                last_pit_index = (last_pit_index + 1) % 14
                if (current_player == 1 and last_pit_index == p2_store_index) or \
                   (current_player == 2 and last_pit_index == p1_store_index):
                    last_pit_index = (last_pit_index + 1) % 14
                board_after_move[last_pit_index] += 1
            
            next_player = 3 - current_player

            # Rule: Go Again
            if (current_player == 1 and last_pit_index == p1_store_index) or \
               (current_player == 2 and last_pit_index == p2_store_index):
                next_player = current_player
            
            # Rule: Capture
            is_own_side = (current_player == 1 and last_pit_index in p1_pits_indices) or \
                          (current_player == 2 and last_pit_index in p2_pits_indices)

            if is_own_side and board_after_move[last_pit_index] == 1:
                opposite_pit_index = 12 - last_pit_index
                if board_after_move[opposite_pit_index] > 0:
                    capture_store = p1_store_index if current_player == 1 else p2_store_index
                    captured_stones = board_after_move[opposite_pit_index] + 1
                    board_after_move[capture_store] += captured_stones
                    board_after_move[opposite_pit_index] = 0
                    board_after_move[last_pit_index] = 0

            new_state = (tuple(board_after_move), next_player)
            if new_state not in visited:
                visited.add(new_state)
                queue.append(new_state)

    # ------------------
    # 4. CONCLUSION
    # ------------------
    print(f"The simulation is complete. The set of all possible score differences is: {sorted(list(final_score_differences))}\n")
    print("Step 3: Final Conclusion\n")
    print("Comparing the possible outcomes to the answer choices:")
    
    choices = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
    unobtainable_choices = []
    
    for choice_letter, diff_value in sorted(choices.items(), key=lambda item: item[1]):
        if diff_value not in final_score_differences:
            unobtainable_choices.append(choice_letter)
            print(f"A score difference of {diff_value} (Choice {choice_letter}) is NOT obtainable.")
        else:
            print(f"A score difference of {diff_value} (Choice {choice_letter}) IS obtainable.")

    print(f"\nThere are {len(unobtainable_choices)} unobtainable score differences in the list: {', '.join([str(choices[c]) for c in unobtainable_choices])}.")
    print("Since more than one of the listed score differences is unobtainable, the correct answer is G.")


# Run the solver
solve_mancala_puzzle()
<<<G>>>