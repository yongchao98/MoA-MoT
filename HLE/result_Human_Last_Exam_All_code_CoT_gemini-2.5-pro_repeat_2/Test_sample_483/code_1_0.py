import functools

# --- Game Mechanics ---

PLAYER_PIECES = {
    1: {"K1", "N1", "R1"},
    2: {"R2", "N2", "K2"}
}

def get_piece_positions(board):
    """Returns a dictionary mapping piece names to their positions."""
    positions = {}
    for i, piece in enumerate(board):
        if piece != ' ':
            positions[piece] = i
    return positions

def is_king_attacked(board, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    p_map = get_piece_positions(board)
    king_key = f'K{player}'
    rook_key = f'R{3 - player}'

    if king_key not in p_map or rook_key not in p_map:
        return False

    k_pos = p_map[king_key]
    r_pos = p_map[rook_key]

    start = min(k_pos, r_pos) + 1
    end = max(k_pos, r_pos)
    
    # Check for any blocking pieces between the king and the rook.
    for i in range(start, end):
        if board[i] != ' ':
            return False  # Path is blocked

    return True

def generate_moves(board, player):
    """Generates all legal next board states for the given player."""
    p_map = get_piece_positions(board)
    my_piece_keys = PLAYER_PIECES[player]
    
    for piece, pos in p_map.items():
        if piece not in my_piece_keys:
            continue

        piece_type = piece[0]
        
        potential_dests = []
        if piece_type == 'K':
            potential_dests.extend([pos - 1, pos + 1])
        elif piece_type == 'N':
            potential_dests.extend([pos - 2, pos + 2])
        elif piece_type == 'R':
            # Move left
            for d in range(pos - 1, -1, -1):
                potential_dests.append(d)
                if board[d] != ' ':
                    break
            # Move right
            for d in range(pos + 1, 8):
                potential_dests.append(d)
                if board[d] != ' ':
                    break
        
        for dest in potential_dests:
            if not (0 <= dest < 8):
                continue

            dest_piece = board[dest]
            if dest_piece != ' ' and dest_piece in my_piece_keys:
                continue

            new_board_list = list(board)
            new_board_list[dest] = piece
            new_board_list[pos] = ' '
            new_board = tuple(new_board_list)
            
            if not is_king_attacked(new_board, player):
                yield new_board

# --- Minimax Solver with Memoization ---

@functools.lru_cache(maxsize=None)
def solve(board, player_to_move):
    """
    Recursively solves the game state using minimax.
    Returns a tuple: (outcome for P1, number of plies).
    'W': P1 wins, 'L': P1 loses, 'D': Draw.
    """
    p_map = get_piece_positions(board)
    if 'K2' not in p_map:
        return 'W', 0
    if 'K1' not in p_map:
        return 'L', 0

    legal_moves = list(generate_moves(board, player_to_move))
    if not legal_moves:
        return 'D', 0

    outcomes = [solve(next_board, 3 - player_to_move) for next_board in legal_moves]

    if player_to_move == 1:
        # P1's turn: Find the best move.
        # 1. Prefer winning moves, and pick the fastest win.
        wins = [ply for res, ply in outcomes if res == 'W']
        if wins:
            return 'W', 1 + min(wins)
        
        # 2. If no win, prefer a draw.
        draws = [ply for res, ply in outcomes if res == 'D']
        if draws:
            return 'D', 0
        
        # 3. If must lose, pick the move that stalls the longest.
        losses = [ply for res, ply in outcomes if res == 'L']
        return 'L', 1 + max(losses)

    else:  # player_to_move == 2
        # P2's turn: Find the best move.
        # 1. Prefer moves that lead to a loss for P1, and pick the fastest.
        losses = [ply for res, ply in outcomes if res == 'L']
        if losses:
            return 'L', 1 + min(losses)

        # 2. If no win for P2, prefer a draw.
        draws = [ply for res, ply in outcomes if res == 'D']
        if draws:
            return 'D', 0
            
        # 3. If P2 must lose, stall for as long as possible.
        wins = [ply for res, ply in outcomes if res == 'W']
        return 'W', 1 + max(wins)

def main():
    """
    Sets up the initial board and solves the game.
    """
    initial_board = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
    
    # Run the solver starting with Player 1's turn
    result, plies = solve(initial_board, 1)

    if result == 'W':
        # A "turn" or move for one player corresponds to 1 ply if it's the first move,
        # or subsequent odd-numbered plies.
        # e.g., plies=1 (P1 wins) -> 1 move
        # e.g., plies=3 (P1, P2, P1 wins) -> 2 moves for P1
        p1_moves = (plies + 1) // 2
        print(p1_moves)
    elif result == 'L':
        print("Player 1 cannot force a win.")
    else: # result == 'D'
        print("The game is a forced draw.")

if __name__ == '__main__':
    main()