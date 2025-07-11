import math

# Memoization cache to store results of solved game states
memo = {}

# --- Game Constants ---
BOARD_SIZE = 8
P1_PIECES = {'K1', 'N1', 'R1'}
P2_PIECES = {'K2', 'N2', 'R2'}

def get_piece_info(piece_name):
    """Gets the owner (1 or 2) and type ('K', 'N', 'R') of a piece."""
    if piece_name == ' ':
        return None, None
    return int(piece_name[1]), piece_name[0]

def find_piece(board, piece_name):
    """Finds the position of a given piece, returning None if not found."""
    try:
        return board.index(piece_name)
    except ValueError:
        return None

def is_king_in_check_by_rook(board, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    king_piece = f'K{player}'
    king_pos = find_piece(board, king_piece)
    if king_pos is None:
        return False

    opponent = 3 - player
    rook_piece = f'R{opponent}'
    rook_pos = find_piece(board, rook_piece)
    if rook_pos is None:
        return False

    start, end = sorted((king_pos, rook_pos))
    for i in range(start + 1, end):
        if board[i] != ' ':
            return False  # Path is blocked
    return True

def apply_move(board, move):
    """Applies a move to the board and returns the new board state as a tuple."""
    from_pos, to_pos = move
    new_board = list(board)
    piece = new_board[from_pos]
    new_board[to_pos] = piece
    new_board[from_pos] = ' '
    return tuple(new_board)

def generate_legal_moves(board, player):
    """Generates all legal moves for the given player."""
    legal_moves = []
    player_pieces = P1_PIECES if player == 1 else P2_PIECES

    for pos, piece in enumerate(board):
        if piece in player_pieces:
            potential_dests = []
            _, piece_type = get_piece_info(piece)

            if piece_type == 'K':
                potential_dests.extend([pos - 1, pos + 1])
            elif piece_type == 'N':
                potential_dests.extend([pos - 2, pos + 2])
            elif piece_type == 'R':
                # Move left
                for i in range(pos - 1, -1, -1):
                    potential_dests.append(i)
                    if board[i] != ' ': break
                # Move right
                for i in range(pos + 1, BOARD_SIZE):
                    potential_dests.append(i)
                    if board[i] != ' ': break

            for dest in potential_dests:
                if 0 <= dest < BOARD_SIZE:
                    dest_owner, _ = get_piece_info(board[dest])
                    if dest_owner != player:
                        temp_board = apply_move(board, (pos, dest))
                        if not is_king_in_check_by_rook(temp_board, player):
                            legal_moves.append((pos, dest))
    return legal_moves

def solve(board, player):
    """
    Recursively solves the game state using minimax with memoization.
    Returns:
        - d > 0: current player can force a win in d plies.
        - float('inf'): game is a draw from this state.
        - d < 0: current player will be forced to lose in abs(d) plies.
    """
    state = (board, player)
    if state in memo:
        return memo[state]

    opponent = 3 - player
    moves = generate_legal_moves(board, player)

    if not moves:
        memo[state] = float('inf')  # Stalemate -> Draw
        return float('inf')

    my_win_plies = []
    opponent_win_plies = []
    is_draw_possible = False

    for move in moves:
        new_board = apply_move(board, move)

        if find_piece(new_board, f'K{opponent}') is None:
            my_win_plies.append(1)
            continue

        res = solve(new_board, opponent)
        
        if res == float('inf'):
            is_draw_possible = True
        elif res < 0:  # Opponent will lose, so I win.
            my_win_plies.append(abs(res) + 1)
        else:  # Opponent will win, so I lose.
            opponent_win_plies.append(res + 1)

    if my_win_plies:
        # If I can win, choose the fastest win.
        result = min(my_win_plies)
    elif is_draw_possible:
        # If I can't win but can draw, do that.
        result = float('inf')
    else:
        # If I must lose, choose the slowest loss.
        result = -max(opponent_win_plies)

    memo[state] = result
    return result

if __name__ == '__main__':
    initial_board = tuple(['K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2'])
    
    # Solve for Player 1 from the start
    plies_to_win = solve(initial_board, 1)

    # Convert plies (half-moves) to Player 1's turns
    # Turns = ceil(plies / 2)
    turns = (int(plies_to_win) + 1) // 2
    
    print(turns)
    print("<<<5>>>")