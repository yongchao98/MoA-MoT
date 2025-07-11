import sys
from functools import lru_cache

# Set a higher recursion limit for deep searches, just in case.
sys.setrecursionlimit(3000)

# Define piece sets and a canonical sorted list for hashing
PIECES_P1 = ['K1', 'N1', 'R1']
PIECES_P2 = ['K2', 'N2', 'R2']
ALL_PIECES_SORTED = sorted(PIECES_P1 + PIECES_P2)
WIN_FOR_P1 = 1
DRAW = 0
WIN_FOR_P2 = -1

def board_to_tuple(board):
    """Converts a board dictionary to a canonical, hashable tuple."""
    return tuple((p, board.get(p)) for p in ALL_PIECES_SORTED)

def is_king_in_check(board, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    king_piece = 'K1' if player == 'P1' else 'K2'
    opp_rook_piece = 'R2' if player == 'P1' else 'R1'

    king_pos = board.get(king_piece)
    rook_pos = board.get(opp_rook_piece)

    if king_pos is None or rook_pos is None:
        return False

    occupied = {pos for pos in board.values() if pos is not None}
    
    start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)
    for pos in range(start + 1, end):
        if pos in occupied:
            return False  # Path is blocked
    return True # Path is clear, king is in check

def get_legal_moves(board, player):
    """Generates all legal moves for a given player from a board state."""
    legal_moves = []
    player_pieces = PIECES_P1 if player == 'P1' else PIECES_P2
    opponent_pieces = PIECES_P2 if player == 'P1' else PIECES_P1

    friendly_pos = {board.get(p) for p in player_pieces if board.get(p) is not None}
    opponent_pos = {board.get(p) for p in opponent_pieces if board.get(p) is not None}

    for piece in player_pieces:
        pos = board.get(piece)
        if pos is None:
            continue

        piece_type = piece[0]
        potential_dests = []

        if piece_type == 'K':
            potential_dests.extend([pos - 1, pos + 1])
        elif piece_type == 'N':
            potential_dests.extend([pos - 2, pos + 2])
        elif piece_type == 'R':
            for i in range(pos - 1, -1, -1): # Left
                potential_dests.append(i)
                if i in friendly_pos or i in opponent_pos: break
            for i in range(pos + 1, 8): # Right
                potential_dests.append(i)
                if i in friendly_pos or i in opponent_pos: break
        
        for dest in potential_dests:
            if not (0 <= dest <= 7) or dest in friendly_pos:
                continue

            new_board = board.copy()
            new_board[piece] = dest
            
            captured_piece = None
            for opp_piece, opp_pos in board.items():
                if opp_pos == dest and opp_piece in opponent_pieces:
                    captured_piece = opp_piece
                    break
            if captured_piece:
                new_board[captured_piece] = None
            
            if not is_king_in_check(new_board, player):
                legal_moves.append(new_board)
    return legal_moves

@lru_cache(maxsize=None)
def solve(board_tuple, player):
    """
    Recursively solves the game state using minimax.
    Returns a tuple: (outcome, plies).
    outcome: 1 for P1 win, 0 for draw, -1 for P2 win.
    plies: number of moves to reach the outcome.
    """
    board = dict(board_tuple)
    moves = get_legal_moves(board, player)

    if not moves:
        # Check for stalemate: no legal moves and opponent's king is alive
        opp_king = 'K2' if player == 'P1' else 'K1'
        if board.get(opp_king) is not None:
            return (DRAW, 0)
        # Should not happen, as king capture would be a move result
        else: # The current player has already lost.
             return (WIN_FOR_P2 if player == 'P1' else WIN_FOR_P1, 0)


    # Initialize best result from the current player's perspective
    if player == 'P1': # Maximizing player
        best_outcome = -2 # Worse than loss
        best_plies = 0
    else: # Player 2, Minimizing player
        best_outcome = 2 # Worse than win
        best_plies = 0

    for next_board in moves:
        # Check for immediate win by capturing the king
        if player == 'P1' and next_board.get('K2') is None:
            return (WIN_FOR_P1, 1)
        if player == 'P2' and next_board.get('K1') is None:
            return (WIN_FOR_P2, 1)
        
        # Recurse to the next state
        next_player = 'P2' if player == 'P1' else 'P1'
        outcome, plies = solve(board_to_tuple(next_board), next_player)
        current_plies = plies + 1

        if player == 'P1':
            if outcome > best_outcome:
                best_outcome = outcome
                best_plies = current_plies
            elif outcome == best_outcome:
                if outcome == WIN_FOR_P1: # If winning, find shortest path
                    best_plies = min(best_plies, current_plies)
                else: # If drawing or losing, find longest path
                    best_plies = max(best_plies, current_plies)
        else: # Player 'P2'
            if outcome < best_outcome:
                best_outcome = outcome
                best_plies = current_plies
            elif outcome == best_outcome:
                if outcome == WIN_FOR_P1: # If P1 wins, P2 stalls (longest path)
                    best_plies = max(best_plies, current_plies)
                else: # If P2 wins or draws, find longest path to seal it
                    best_plies = max(best_plies, current_plies)

    return (best_outcome, best_plies)


def find_forced_win_turns():
    """Sets up the initial board and starts the solver."""
    initial_board = {
        'K1': 0, 'N1': 1, 'R1': 2,
        'K2': 7, 'N2': 6, 'R2': 5,
    }
    initial_board_tuple = board_to_tuple(initial_board)
    
    outcome, plies = solve(initial_board_tuple, 'P1')

    if outcome == WIN_FOR_P1:
        print(plies)
    else:
        print("Player 1 cannot force a win from the initial position.")


if __name__ == "__main__":
    find_forced_win_turns()