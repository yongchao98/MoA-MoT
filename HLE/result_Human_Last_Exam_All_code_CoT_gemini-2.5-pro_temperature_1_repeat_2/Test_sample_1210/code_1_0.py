import collections

def solve_mancala():
    """
    Simulates a Mancala game from a specific state to find all possible score differences.
    """
    final_diffs = set()

    # Board representation:
    # A tuple of 14 integers.
    # board[0-5]: Player 1's pits
    # board[6]: Player 1's store
    # board[7-12]: Player 2's pits
    # board[13]: Player 2's store
    
    # Initial state from the problem description
    p1_pits = (0, 2, 0, 0, 2, 0)
    p1_store = 22
    # P2 pits are given left-to-right from P2's perspective,
    # which corresponds to indices 7 through 12 in our model.
    p2_pits = (1, 0, 0, 0, 0, 0)
    p2_store = 21
    
    initial_board = p1_pits + (p1_store,) + p2_pits + (p2_store,)
    
    # state is (board_tuple, turn_player_number)
    q = collections.deque([(initial_board, 1)])
    visited = {(initial_board, 1)}

    while q:
        board_t, turn = q.popleft()
        
        # Check for game over condition
        p1_pits_empty = sum(board_t[0:6]) == 0
        p2_pits_empty = sum(board_t[7:13]) == 0
        
        if p1_pits_empty or p2_pits_empty:
            final_board = list(board_t)
            # Collect remaining stones
            if p1_pits_empty:
                final_board[13] += sum(final_board[7:13])
                for i in range(7, 13): final_board[i] = 0
            if p2_pits_empty:
                final_board[6] += sum(final_board[0:6])
                for i in range(0, 6): final_board[i] = 0
            
            diff = abs(final_board[6] - final_board[13])
            final_diffs.add(diff)
            continue
        
        if turn == 1:
            pit_indices = [i for i in range(6) if board_t[i] > 0]
            if not pit_indices: # Player must skip turn if no moves
                if (board_t, 2) not in visited:
                    visited.add((board_t, 2))
                    q.append((board_t, 2))
                continue
            
            for pit_idx in pit_indices:
                board = list(board_t)
                stones = board[pit_idx]
                board[pit_idx] = 0
                
                current_idx = pit_idx
                for _ in range(stones):
                    current_idx = (current_idx + 1) % 14
                    if current_idx == 13: current_idx = 0 # Skip opponent's store
                    board[current_idx] += 1
                
                # Turn logic: Go Again or Capture
                if current_idx == 6: # Go Again
                    next_turn = 1
                elif 0 <= current_idx <= 5 and board[current_idx] == 1: # Possible Capture
                    opposite_idx = 12 - current_idx
                    if board[opposite_idx] > 0:
                        board[6] += 1 + board[opposite_idx]
                        board[current_idx] = 0
                        board[opposite_idx] = 0
                    next_turn = 2
                else: # Normal turn end
                    next_turn = 2
                
                if (tuple(board), next_turn) not in visited:
                    visited.add((tuple(board), next_turn))
                    q.append((tuple(board), next_turn))

        else: # turn == 2
            pit_indices = [i for i in range(7, 13) if board_t[i] > 0]
            if not pit_indices: # Player must skip turn
                if (board_t, 1) not in visited:
                    visited.add((board_t, 1))
                    q.append((board_t, 1))
                continue
            
            for pit_idx in pit_indices:
                board = list(board_t)
                stones = board[pit_idx]
                board[pit_idx] = 0
                
                current_idx = pit_idx
                for _ in range(stones):
                    current_idx = (current_idx + 1) % 14
                    if current_idx == 6: current_idx = 7 # Skip opponent's store
                    board[current_idx] += 1

                if current_idx == 13: # Go Again
                    next_turn = 2
                elif 7 <= current_idx <= 12 and board[current_idx] == 1: # Possible Capture
                    opposite_idx = 12 - current_idx
                    if board[opposite_idx] > 0:
                        board[13] += 1 + board[opposite_idx]
                        board[current_idx] = 0
                        board[opposite_idx] = 0
                    next_turn = 1
                else: # Normal turn end
                    next_turn = 1
                
                if (tuple(board), next_turn) not in visited:
                    visited.add((tuple(board), next_turn))
                    q.append((tuple(board), next_turn))
    
    possible_differences = sorted(list(final_diffs))
    choices = {"Zero": 0, "One": 1, "Two": 2, "Three": 3, "Four": 4, "Five": 5}
    
    print("--- Mancala Game Simulation ---")
    print(f"Initial P1 Pits: {p1_pits}, Store: {p1_store}")
    print(f"Initial P2 Pits: {p2_pits}, Store: {p2_store}")
    print("\nSimulating all game paths...")
    print(f"\nReachable final score differences: {possible_differences}")
    print("\nChecking against answer choices:")
    for name, value in choices.items():
        if value in possible_differences:
            print(f"- {name} ({value}): Possible")
        else:
            print(f"- {name} ({value}): Not Possible")
            
solve_mancala()