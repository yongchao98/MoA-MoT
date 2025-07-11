from collections import deque

def solve_mancala():
    """
    Simulates all possible game outcomes from the given state to find all
    possible score differences.
    """
    # Initial state from the problem description
    # p1_pits[i] is opposite p2_pits[5-i]
    initial_p1_pits = (0, 2, 0, 0, 2, 0)
    initial_p1_store = 22
    initial_p2_pits = (1, 0, 0, 0, 0, 0)
    initial_p2_store = 21
    
    # Queue for breadth-first search of game states
    # State: (p1_pits, p2_pits, p1_store, p2_store, player_turn)
    q = deque([(initial_p1_pits, initial_p2_pits, initial_p1_store, initial_p2_store, 1)])
    
    # Set to keep track of visited states to avoid loops and redundant computations
    visited = set()
    
    # Set to store the unique final score differences found
    final_diffs = set()
    
    print("Simulating all possible game outcomes...\n")

    while q:
        p1_pits, p2_pits, p1_store, p2_store, turn = q.popleft()

        state_key = (p1_pits, p2_pits, p1_store, p2_store, turn)
        if state_key in visited:
            continue
        visited.add(state_key)

        # Check for game over condition
        if sum(p1_pits) == 0 or sum(p2_pits) == 0:
            final_p1_score = p1_store + sum(p1_pits)
            final_p2_score = p2_store + sum(p2_pits)
            
            winner_score = max(final_p1_score, final_p2_score)
            loser_score = min(final_p1_score, final_p2_score)
            
            diff = winner_score - loser_score
            
            if diff not in final_diffs:
                final_diffs.add(diff)
                print(f"Game end found: P1 Score={final_p1_score}, P2 Score={final_p2_score}")
                print(f"Final equation: {winner_score} - {loser_score} = {diff}\n")
            continue

        # Generate next moves for the current player
        if turn == 1:
            possible_moves = [i for i, s in enumerate(p1_pits) if s > 0]
            for move_idx in possible_moves:
                p1p, p2p, p1s, p2s = list(p1_pits), list(p2_pits), p1_store, p2_store
                
                stones = p1p[move_idx]
                p1p[move_idx] = 0
                
                pit = move_idx
                while stones > 0:
                    pit = (pit + 1) % 13  # P1 sows through P1 pits, P1 store, P2 pits (13 spots)
                    if pit < 6: p1p[pit] += 1
                    elif pit == 6: p1s += 1
                    else: p2p[pit - 7] += 1
                    stones -= 1

                # Check Go Again
                if pit == 6:
                    q.append((tuple(p1p), tuple(p2p), p1s, p2s, 1))
                    continue
                
                # Check Capture
                if pit < 6 and p1p[pit] == 1:
                    opposite_pit = 5 - pit
                    if p2p[opposite_pit] > 0:
                        p1s += p2p[opposite_pit] + 1
                        p2p[opposite_pit] = 0
                        p1p[pit] = 0

                q.append((tuple(p1p), tuple(p2p), p1s, p2s, 2))
        else:  # turn == 2
            possible_moves = [i for i, s in enumerate(p2_pits) if s > 0]
            for move_idx in possible_moves:
                p1p, p2p, p1s, p2s = list(p1_pits), list(p2_pits), p1_store, p2_store
                
                stones = p2p[move_idx]
                p2p[move_idx] = 0
                
                # Sowing for P2 uses a different path on the board representation
                board_idx = move_idx + 7
                while stones > 0:
                    board_idx = (board_idx + 1) % 14
                    if board_idx == 6: continue # Skip P1's store
                    
                    if board_idx < 6: p1p[board_idx] += 1
                    elif 7 <= board_idx < 13: p2p[board_idx - 7] += 1
                    elif board_idx == 13: p2s += 1
                    stones -= 1
                
                # Check Go Again
                if board_idx == 13:
                    q.append((tuple(p1p), tuple(p2p), p1s, p2s, 2))
                    continue
                
                # Check Capture
                if 7 <= board_idx < 13 and p2p[board_idx - 7] == 1:
                    opposite_pit = 5 - (board_idx - 7)
                    if p1p[opposite_pit] > 0:
                        p2s += p1p[opposite_pit] + 1
                        p1p[opposite_pit] = 0
                        p2p[board_idx - 7] = 0

                q.append((tuple(p1p), tuple(p2p), p1s, p2s, 1))

    print("-----------------------------------------")
    print(f"Possible score differences are: {sorted(list(final_diffs))}")
    
    choices = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
    impossible_choices = []
    for choice, value in choices.items():
        if value not in final_diffs:
            impossible_choices.append(f"{choice}. {value}")
    
    print(f"The score differences from the list that are NOT possible are: {impossible_choices}")
    print("Based on the analysis, all odd score differences are impossible.")

solve_mancala()