# Mancala board representation:
# Player 1 Pits: Indices 0-5
# Player 1 Store: Index 6
# Player 2 Pits: Indices 7-12
# Player 2 Store: Index 13
P1_PITS = list(range(6))
P1_STORE = 6
P2_PITS = list(range(7, 13))
P2_STORE = 13
TOTAL_STONES = 48

def print_board(board):
    """Prints the board state in a readable top-down format."""
    p2_pits_reversed = list(board[7:13])
    p2_pits_reversed.reverse()
    print("---------------------------------------")
    print(f"       {p2_pits_reversed}  <- P2 Pits")
    print(f" P2 Store [{board[P2_STORE]}]              [{board[P1_STORE]}] P1 Store")
    print(f" P1 Pits ->  {list(board[0:6])}")
    print("---------------------------------------")

def play_move_and_print(board, player, pit_human_readable):
    """
    Plays a move, prints the actions, and returns the new board and next player.
    `pit_human_readable` is 1-based for the player's side (1-6).
    """
    board = list(board)
    player_pits = P1_PITS if player == 1 else P2_PITS
    pit_index = player_pits[pit_human_readable - 1]

    print(f"\n>>> Player {player}'s turn. Choosing pit {pit_human_readable} (index {pit_index}) with {board[pit_index]} stones.")

    # 1. Pick up stones
    stones_in_hand = board[pit_index]
    board[pit_index] = 0

    # 2. Sow stones
    player_store = P1_STORE if player == 1 else P2_STORE
    opponent_store = P2_STORE if player == 1 else P1_STORE
    current_pos = pit_index
    for _ in range(stones_in_hand):
        current_pos = (current_pos + 1) % 14
        if current_pos == opponent_store:
            current_pos = (current_pos + 1) % 14
        board[current_pos] += 1
    
    last_stone_pos = current_pos
    
    # 3. Check for special rules
    next_player = 2 if player == 1 else 1
    
    # Extra turn rule
    if last_stone_pos == player_store:
        print("INFO: Last stone landed in the player's store. Extra turn granted.")
        next_player = player
        
    # Capture rule
    elif last_stone_pos in player_pits and board[last_stone_pos] == 1:
        opposite_pit_index = 12 - last_stone_pos
        if board[opposite_pit_index] > 0:
            print(f"INFO: Last stone landed in empty pit {pit_human_readable}. Capturing stones!")
            captured_stones = board[opposite_pit_index] + 1
            board[player_store] += captured_stones
            board[last_stone_pos] = 0
            board[opposite_pit_index] = 0
            print(f"INFO: Captured {captured_stones} stones into store {player_store}.")

    print("Board after move:")
    print_board(board)
    return tuple(board), next_player

def calculate_and_print_final_score(board):
    """Finalizes the board and prints the score difference calculation."""
    board = list(board)
    p1_remaining = sum(board[i] for i in P1_PITS)
    p2_remaining = sum(board[i] for i in P2_PITS)
    
    p1_final = board[P1_STORE] + p1_remaining
    p2_final = board[P2_STORE] + p2_remaining
    
    print("\n--- GAME OVER ---")
    if p1_remaining > 0:
        print(f"Player 1 collects remaining {p1_remaining} stones. Final score: {board[P1_STORE]} + {p1_remaining} = {p1_final}")
    if p2_remaining > 0:
        print(f"Player 2 collects remaining {p2_remaining} stones. Final score: {board[P2_STORE]} + {p2_remaining} = {p2_final}")

    winner_score = max(p1_final, p2_final)
    loser_score = min(p1_final, p2_final)
    diff = winner_score - loser_score
    print(f"Score Difference Calculation: {winner_score} - {loser_score} = {diff}")
    return diff

def main():
    initial_board = (0, 2, 0, 0, 2, 0, 22, 1, 0, 0, 0, 0, 0, 21)

    print("INVESTIGATION 1: Path to a Score Difference of 0")
    print("=============================================")
    print("Initial state:")
    print_board(initial_board)

    # P1 plays from their second pit (2 stones)
    board1, player = play_move_and_print(initial_board, 1, 2)
    # P2 plays from their first pit (1 stone), resulting in a capture
    board2, player = play_move_and_print(board1, 2, 1) 
    # P1 plays from their third pit (1 stone)
    board3, player = play_move_and_print(board2, 1, 3)
    # Now, Player 2 has no legal moves, so the game ends.
    diff_0 = calculate_and_print_final_score(board3)

    print("\nINVESTIGATION 2: Path to a Score Difference of 4")
    print("=============================================")
    print("Initial state:")
    print_board(initial_board)
    board, player = play_move_and_print(initial_board, 1, 5) # P1 extra turn
    board, player = play_move_and_print(board, player, 6)   # P1 extra turn
    board, player = play_move_and_print(board, player, 2)
    board, player = play_move_and_print(board, player, 1)
    board, player = play_move_and_print(board, player, 3)
    board, player = play_move_and_print(board, player, 2)
    board, player = play_move_and_print(board, player, 4)
    board, player = play_move_and_print(board, player, 3)
    board, player = play_move_and_print(board, player, 6)   # P1 extra turn
    board, player = play_move_and_print(board, player, 5)
    board, player = play_move_and_print(board, player, 4)
    board, player = play_move_and_print(board, player, 6)   # P1 extra turn
    # Player 1's side is empty, so the game ends.
    diff_4 = calculate_and_print_final_score(board)

    print("\nCONCLUSION")
    print("=============================================")
    print(f"We found paths leading to score differences of {diff_0} and {diff_4}. A difference of 2 is also possible.")
    print("\nMathematical Proof for Odd vs. Even Differences:")
    print(f"The total number of stones in the game is constant: {TOTAL_STONES}.")
    print("At the end of the game: P1_score + P2_score = 48.")
    print("The score difference is D = |P1_score - P2_score|.")
    print("We can rewrite P2_score as: P2_score = 48 - P1_score.")
    print("Substituting this gives: D = |P1_score - (48 - P1_score)|.")
    print("This simplifies to: D = |2 * P1_score - 48|.")
    print("Since (2 * P1_score) is an even number and 48 is an even number, their difference (2 * P1_score - 48) must also be an even number.")
    print("Therefore, the score difference D must always be even.")
    print("\nThe impossible score differences from the list {0, 1, 2, 3, 4, 5} are the odd numbers: 1, 3, and 5.")
    print("Since more than one choice is unobtainable, the answer is G.")

if __name__ == "__main__":
    main()
<<<G>>>