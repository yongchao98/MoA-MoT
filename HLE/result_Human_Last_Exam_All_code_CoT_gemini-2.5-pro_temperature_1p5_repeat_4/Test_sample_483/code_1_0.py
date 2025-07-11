import sys
import collections

# Increase recursion limit for potentially deep game trees.
sys.setrecursionlimit(2000)

# --- Game Configuration ---
P1_PIECES = {'K1', 'N1', 'R1'}
P2_PIECES = {'K2', 'N2', 'R2'}
PIECE_OWNER = {p: 1 for p in P1_PIECES} | {p: 2 for p in P2_PIECES}
memo = {} # Memoization table for storing state outcomes.

def get_piece_positions(board):
    """Returns a dictionary mapping piece names to their board indices."""
    return {piece: i for i, piece in enumerate(board) if piece != ' '}

def is_king_in_check(board, player):
    """Checks if the specified player's King is under attack by the opponent's Rook."""
    positions = get_piece_positions(board)
    king_piece = f'K{player}'
    opponent_rook_piece = f'R{3 - player}'

    if king_piece not in positions or opponent_rook_piece not in positions:
        return False

    k_pos, r_pos = positions[king_piece], positions[opponent_rook_piece]
    start, end = (k_pos + 1, r_pos) if k_pos < r_pos else (r_pos + 1, k_pos)
    
    # Check if the path between them is clear of any other pieces.
    return all(board[i] == ' ' for i in range(start, end))

def get_legal_moves(board, player):
    """Generates a list of all legal moves (from_pos, to_pos) for the given player."""
    legal_moves = []
    positions = get_piece_positions(board)
    
    for piece, pos in positions.items():
        if PIECE_OWNER.get(piece) != player:
            continue

        potential_dests = []
        piece_type = piece[0]
        if piece_type == 'K': potential_dests.extend([pos - 1, pos + 1])
        elif piece_type == 'N': potential_dests.extend([pos - 2, pos + 2])
        elif piece_type == 'R':
            for direction in [-1, 1]:
                for i in range(pos + direction, 8 if direction == 1 else -1, direction):
                    potential_dests.append(i)
                    if board[i] != ' ': break
        
        for dest in potential_dests:
            if not (0 <= dest <= 7) or (board[dest] != ' ' and PIECE_OWNER.get(board[dest]) == player):
                continue
            
            next_board = list(board)
            next_board[dest], next_board[pos] = piece, ' '
            if not is_king_in_check(tuple(next_board), player):
                legal_moves.append((pos, dest))
    return legal_moves

def find_outcome(board, player):
    """
    Recursively determines the game outcome using minimax with memoization.
    Returns (outcome, half_moves) for the current player.
    """
    state = (board, player)
    if state in memo: return memo[state]

    positions = get_piece_positions(board)
    if f'K{player}' not in positions: return ('LOSS', 0)
    if f'K{3-player}' not in positions: return ('WIN', 0)

    moves = get_legal_moves(board, player)
    if not moves:
        return ('LOSS', 0) if is_king_in_check(board, player) else ('DRAW', float('inf'))

    can_force_win, can_force_draw = False, False
    min_win_turns, max_loss_turns = float('inf'), 0

    for pos, dest in moves:
        next_board_list = list(board)
        next_board_list[dest], next_board_list[pos] = board[pos], ' '
        outcome, turns = find_outcome(tuple(next_board_list), 3 - player)
        
        if outcome == 'LOSS':
            can_force_win = True
            min_win_turns = min(min_win_turns, turns + 1)
        elif outcome == 'DRAW':
            can_force_draw = True
        else: # WIN
            max_loss_turns = max(max_loss_turns, turns + 1)

    if can_force_win: result = ('WIN', min_win_turns)
    elif can_force_draw: result = ('DRAW', float('inf'))
    else: result = ('LOSS', max_loss_turns)
        
    memo[state] = result
    return result

def solve_and_print():
    """Solves the game and prints the result in the required format."""
    initial_board = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
    outcome, half_moves = find_outcome(initial_board, 1)

    if outcome == 'WIN':
        player1_turns = (half_moves + 1) // 2
        print(f"Player 1 can force a win in a minimum of {player1_turns} moves.")
        print("The number of Player 1's moves is calculated from the total half-moves (plies):")
        print(f"({half_moves} + 1) // 2 = {player1_turns}")
    elif outcome == 'DRAW':
        print("The game results in a draw with optimal play.")
    else:
        print("Player 1 cannot force a win against an optimal opponent.")

if __name__ == "__main__":
    solve_and_print()