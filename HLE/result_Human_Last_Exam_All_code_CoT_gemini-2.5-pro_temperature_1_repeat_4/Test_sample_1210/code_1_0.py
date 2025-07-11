from collections import deque

def solve_mancala():
    """
    Simulates all possible game outcomes from the given Mancala state
    and determines which score differences are possible.
    """
    # Board layout: p1_pits[0-5], p1_store[6], p2_pits[7-12], p2_store[13]
    initial_board = (0, 2, 0, 0, 2, 0, 22, 1, 0, 0, 0, 0, 0, 21)
    
    # Use a queue for Breadth-First Search of game states
    # Item format: (board_tuple, player_turn) where player 1 is 0, player 2 is 1
    queue = deque([(initial_board, 0)])
    
    # Keep track of visited states to avoid cycles and redundant computations
    visited = set()
    
    # Store final score outcomes found
    final_outcomes = set()

    while queue:
        current_board_tuple, player_turn = queue.popleft()

        if current_board_tuple in visited:
            continue
        visited.add(current_board_tuple)

        board = list(current_board_tuple)

        # Determine pit ranges for the current player
        if player_turn == 0: # Player 1
            pit_indices = range(0, 6)
            player_store_idx = 6
            opponent_store_idx = 13
        else: # Player 2
            pit_indices = range(7, 13)
            player_store_idx = 13
            opponent_store_idx = 6

        # Find valid moves (pits with stones)
        possible_moves = [i for i in pit_indices if board[i] > 0]

        for move_pit_idx in possible_moves:
            temp_board = list(board)
            
            stones_to_sow = temp_board[move_pit_idx]
            temp_board[move_pit_idx] = 0
            
            current_pos = move_pit_idx
            while stones_to_sow > 0:
                current_pos = (current_pos + 1) % 14
                if current_pos == opponent_store_idx:
                    continue # Skip opponent's store
                
                temp_board[current_pos] += 1
                stones_to_sow -= 1
            
            last_stone_pos = current_pos
            
            # Capture rule
            # Check if last stone landed in player's own empty pit
            is_own_pit = last_stone_pos in pit_indices
            was_empty = temp_board[last_stone_pos] == 1
            if is_own_pit and was_empty:
                opposite_pit_idx = 12 - last_stone_pos
                if temp_board[opposite_pit_idx] > 0:
                    captured_stones = temp_board[opposite_pit_idx] + temp_board[last_stone_pos]
                    temp_board[player_store_idx] += captured_stones
                    temp_board[opposite_pit_idx] = 0
                    temp_board[last_stone_pos] = 0
            
            # Check for game end
            p1_pits_empty = all(s == 0 for s in temp_board[0:6])
            p2_pits_empty = all(s == 0 for s in temp_board[7:13])
            
            if p1_pits_empty or p2_pits_empty:
                # Sweep remaining stones
                p1_rem = sum(temp_board[0:6])
                p2_rem = sum(temp_board[7:13])
                temp_board[6] += p1_rem
                temp_board[13] += p2_rem
                for i in range(0,6): temp_board[i] = 0
                for i in range(7,13): temp_board[i] = 0
                
                final_outcomes.add((temp_board[6], temp_board[13]))
                continue

            # "Go Again" rule
            if last_stone_pos == player_store_idx:
                next_player = player_turn
            else:
                next_player = 1 - player_turn # Switch player
            
            queue.append((tuple(temp_board), next_player))

    print("All possible final scores (P1, P2) and differences:")
    possible_differences = set()
    for s1, s2 in sorted(list(final_outcomes)):
        winner = max(s1, s2)
        loser = min(s1, s2)
        diff = winner - loser
        possible_differences.add(diff)
        print(f"P1: {s1}, P2: {s2} -> Score difference: {winner} - {loser} = {diff}")
        
    print("\n--- Analysis ---")
    print(f"The set of all possible score differences is: {sorted(list(possible_differences))}")
    
    choices = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
    unobtainable = []
    for choice, value in choices.items():
        if value not in possible_differences:
            unobtainable.append(value)
    
    print(f"The unobtainable score differences from the list are: {unobtainable}")
    
    if len(unobtainable) > 1:
        print("More than one of the listed score differences is unobtainable.")
    elif len(unobtainable) == 1:
        print(f"The only unobtainable score difference is {unobtainable[0]}.")
    else:
        print("All listed score differences are obtainable.")


solve_mancala()
<<<G>>>