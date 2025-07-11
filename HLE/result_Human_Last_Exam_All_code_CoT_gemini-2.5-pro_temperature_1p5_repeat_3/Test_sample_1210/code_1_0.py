import collections

# Using a tuple for board state so it can be used as a dictionary key for memoization
# Format: (p1_pits_tuple, p1_store, p2_pits_tuple, p2_store)
initial_state = ([0, 2, 0, 0, 2, 0], 22, [1, 0, 0, 0, 0, 0], 21)
player_to_move = 1

memo = {}
final_score_diffs = set()

def get_possible_moves(pits, player):
    """Returns a list of valid pit indices to play from."""
    return [i for i, stones in enumerate(pits) if stones > 0]

def search(p1_pits, p1_store, p2_pits, p2_store, player):
    """
    Recursively explores all game paths from the current state.
    """
    state_key = (tuple(p1_pits), p1_store, tuple(p2_pits), p2_store, player)
    if state_key in memo:
        return
    
    # Check for game over condition
    if sum(p1_pits) == 0 or sum(p2_pits) == 0:
        final_p1_score = p1_store + sum(p1_pits)
        final_p2_score = p2_store + sum(p2_pits)
        diff = abs(final_p1_score - final_p2_score)
        if diff not in final_score_diffs:
            print(f"Game ended. Player 1: {final_p1_score}, Player 2: {final_p2_score}. Difference: {diff}")
            final_score_diffs.add(diff)
        return

    memo[state_key] = True

    if player == 1:
        moves = get_possible_moves(p1_pits, player)
        for move_pit in moves:
            # Create a copy of the state for this move
            temp_p1_pits, temp_p1_store = list(p1_pits), p1_store
            temp_p2_pits, temp_p2_store = list(p2_pits), p2_store

            stones = temp_p1_pits[move_pit]
            temp_p1_pits[move_pit] = 0
            
            current_side_pits = temp_p1_pits
            current_pit_idx = move_pit

            for i in range(stones):
                current_pit_idx += 1
                if current_pit_idx == 6: # P1's store
                    temp_p1_store += 1
                    if i == stones - 1: # Last stone
                        search(temp_p1_pits, temp_p1_store, temp_p2_pits, temp_p2_store, 1) # Go again
                        continue
                elif current_pit_idx > 6: # Switch to P2's side
                    current_side_pits = temp_p2_pits
                    current_pit_idx = 0
                    current_side_pits[current_pit_idx] += 1
                else: # Still on P1's side
                    current_side_pits[current_pit_idx] += 1
            
            # If not a "go again" move
            else:
                last_pit_idx = current_pit_idx
                # Capture Rule check
                if current_side_pits == temp_p1_pits and temp_p1_pits[last_pit_idx] == 1:
                    opposite_pit_idx = 5 - last_pit_idx
                    if temp_p2_pits[opposite_pit_idx] > 0:
                        captured_stones = temp_p2_pits[opposite_pit_idx] + temp_p1_pits[last_pit_idx]
                        temp_p2_pits[opposite_pit_idx] = 0
                        temp_p1_pits[last_pit_idx] = 0
                        temp_p1_store += captured_stones
                
                search(temp_p1_pits, temp_p1_store, temp_p2_pits, temp_p2_store, 2)
    else: # Player 2's turn
        moves = get_possible_moves(p2_pits, player)
        for move_pit in moves:
            # Create a copy of the state for this move
            temp_p1_pits, temp_p1_store = list(p1_pits), p1_store
            temp_p2_pits, temp_p2_store = list(p2_pits), p2_store

            stones = temp_p2_pits[move_pit]
            temp_p2_pits[move_pit] = 0
            
            current_side_pits = temp_p2_pits
            current_pit_idx = move_pit

            for i in range(stones):
                current_pit_idx += 1
                if current_pit_idx == 6: # P2's store
                    temp_p2_store += 1
                    if i == stones - 1: # Last stone
                        search(temp_p1_pits, temp_p1_store, temp_p2_pits, temp_p2_store, 2) # Go again
                        continue
                elif current_pit_idx > 6: # Switch to P1's side
                    current_side_pits = temp_p1_pits
                    current_pit_idx = 0
                    current_side_pits[current_pit_idx] += 1
                else: # Still on P2's side
                    current_side_pits[current_pit_idx] += 1
            
            # If not a "go again" move
            else:
                last_pit_idx = current_pit_idx
                # Capture Rule check
                if current_side_pits == temp_p2_pits and temp_p2_pits[last_pit_idx] == 1:
                    opposite_pit_idx = 5 - last_pit_idx
                    if temp_p1_pits[opposite_pit_idx] > 0:
                        captured_stones = temp_p1_pits[opposite_pit_idx] + temp_p2_pits[last_pit_idx]
                        temp_p1_pits[opposite_pit_idx] = 0
                        temp_p2_pits[last_pit_idx] = 0
                        temp_p2_store += captured_stones
                
                search(temp_p1_pits, temp_p1_store, temp_p2_pits, temp_p2_store, 1)

def solve():
    print("Simulating all possible game outcomes from the initial state...")
    search(initial_state[0], initial_state[1], initial_state[2], initial_state[3], player_to_move)

    print("\n----------------------------------------------------")
    print(f"All possible score differences found: {sorted(list(final_score_diffs))}")
    print("----------------------------------------------------")

    choices = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5}
    unobtainable_choices = []

    print("Analyzing the answer choices:")
    for choice, diff in choices.items():
        if diff in final_score_diffs:
            print(f"{choice}. {diff}: This score difference is possible.")
        else:
            print(f"{choice}. {diff}: This score difference is NOT possible.")
            unobtainable_choices.append(choice)
    
    print("\nFinal Conclusion:")
    if len(unobtainable_choices) > 1:
        print("More than one of the listed score differences is unobtainable.")
    else:
        print(f"The only unobtainable score difference is {unobtainable_choices[0]}.")
    
    total_stones = sum(initial_state[0]) + initial_state[1] + sum(initial_state[2]) + initial_state[3]
    print(f"\nMathematical Check: The total number of stones is {total_stones}.")
    print("Since the total is an even number, the difference between final scores must also be even.")
    print("This confirms that any odd difference (1, 3, 5) is impossible.")

solve()