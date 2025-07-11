from collections import deque

def solve_mancala():
    """
    Solves the Mancala problem by exploring all possible game states.
    """
    # Initial game state
    # P1 pits (indices 0-5), P1 store (index 6)
    # P2 pits (indices 7-12), P2 store (index 13)
    initial_board = [0, 2, 0, 0, 2, 0, 22, 1, 0, 0, 0, 0, 0, 21]
    initial_turn = 1  # Player 1's turn

    # Use a deque for Breadth-First Search (BFS)
    # The queue stores tuples of (board_state, current_turn)
    queue = deque([(tuple(initial_board), initial_turn)])
    
    # A set to keep track of visited states to avoid cycles and redundant work
    # A state is defined by the board and whose turn it is.
    visited = { (tuple(initial_board), initial_turn) }
    
    # A set to store the final score differences found
    found_differences = set()

    while queue:
        current_board_tuple, turn = queue.popleft()
        current_board = list(current_board_tuple)

        # Check for game-ending conditions
        p1_pits_empty = sum(current_board[0:6]) == 0
        p2_pits_empty = sum(current_board[7:12]) == 0

        if p1_pits_empty or p2_pits_empty:
            final_board = list(current_board)
            # Player 2 collects remaining stones
            final_board[13] += sum(final_board[7:12])
            for i in range(7, 13):
                final_board[i] = 0
            # Player 1 collects remaining stones
            final_board[6] += sum(final_board[0:6])
            for i in range(0, 6):
                final_board[i] = 0
            
            p1_score = final_board[6]
            p2_score = final_board[13]
            diff = abs(p1_score - p2_score)
            
            if diff not in found_differences:
                found_differences.add(diff)
                winner_score = max(p1_score, p2_score)
                loser_score = min(p1_score, p2_score)
                print(f"Found possible outcome. Scores: {winner_score} vs {loser_score}.")
                print(f"Final Equation: {winner_score} - {loser_score} = {diff}")
            continue

        # Generate next possible states from moves
        if turn == 1:
            pit_indices = range(0, 6)
            next_turn = 2
        else:  # turn == 2
            pit_indices = range(7, 13)
            next_turn = 1
        
        for pit_idx in pit_indices:
            if current_board[pit_idx] == 0:
                continue # Cannot move from an empty pit

            board_after_move = list(current_board)
            stones = board_after_move[pit_idx]
            board_after_move[pit_idx] = 0
            
            # Distribute stones
            current_pos = pit_idx
            for _ in range(stones):
                current_pos = (current_pos + 1) % 14
                # Skip opponent's store
                if turn == 1 and current_pos == 13: # P1 skips P2's store
                    current_pos = 0
                if turn == 2 and current_pos == 6: # P2 skips P1's store
                    current_pos = 7
                board_after_move[current_pos] += 1
            
            last_stone_pos = current_pos
            
            # Capture rule
            is_capture = False
            # P1 Capture
            if turn == 1 and 0 <= last_stone_pos <= 5 and board_after_move[last_stone_pos] == 1:
                opposite_pit_idx = 12 - last_stone_pos
                if board_after_move[opposite_pit_idx] > 0:
                    board_after_move[6] += board_after_move[last_stone_pos] + board_after_move[opposite_pit_idx]
                    board_after_move[last_stone_pos] = 0
                    board_after_move[opposite_pit_idx] = 0
                    is_capture = True
            # P2 Capture
            elif turn == 2 and 7 <= last_stone_pos <= 12 and board_after_move[last_stone_pos] == 1:
                opposite_pit_idx = 12 - last_stone_pos
                if board_after_move[opposite_pit_idx] > 0:
                    board_after_move[13] += board_after_move[last_stone_pos] + board_after_move[opposite_pit_idx]
                    board_after_move[last_stone_pos] = 0
                    board_after_move[opposite_pit_idx] = 0
                    is_capture = True

            # "Go Again" rule
            turn_for_next_state = next_turn
            if not is_capture:
                if (turn == 1 and last_stone_pos == 6) or (turn == 2 and last_stone_pos == 13):
                    turn_for_next_state = turn

            new_state = (tuple(board_after_move), turn_for_next_state)
            if new_state not in visited:
                visited.add(new_state)
                queue.append(new_state)

    print("\n--- Analysis ---")
    print(f"All possible score differences found: {sorted(list(found_differences))}")

    choices = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
    unobtainable = []
    
    print("\nChecking answer choices:")
    for choice, diff_val in choices.items():
        if diff_val not in found_differences:
            unobtainable.append(choice)
            print(f"Score difference of {diff_val} is NOT obtainable.")
    
    if len(unobtainable) > 1:
        print("\nMore than one of the listed score differences is unobtainable.")
    
if __name__ == '__main__':
    solve_mancala()
    print("\n<<<G>>>")