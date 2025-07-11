import sys

# It's a deep recursion problem, increase the limit.
sys.setrecursionlimit(2000)

# Global memoization cache
memo = {}

# Piece definitions
PLAYER_1_PIECES = {'K1', 'N1', 'R1'}
PLAYER_2_PIECES = {'K2', 'N2', 'R2'}

# --- Helper Functions ---

def get_piece_positions(board):
    """Converts board tuple to a dictionary of piece positions."""
    positions = {}
    for i, piece in enumerate(board):
        if piece != ' ':
            positions[piece] = i
    return positions

def apply_move(board, move):
    """Applies a move and returns the new board tuple."""
    from_pos, to_pos = move
    list_board = list(board)
    piece = list_board[from_pos]
    list_board[to_pos] = piece
    list_board[from_pos] = ' '
    return tuple(list_board)

def is_in_check(board, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    positions = get_piece_positions(board)
    player_king = f'K{player}'
    opponent_rook = f'R{3 - player}'

    if player_king not in positions or opponent_rook not in positions:
        return False

    king_pos = positions[player_king]
    rook_pos = positions[opponent_rook]

    start, end = sorted((king_pos, rook_pos))
    for i in range(start + 1, end):
        if board[i] != ' ':
            return False  # Path is blocked

    return True

def get_legal_moves(board, player):
    """Generates all legal moves for a given player, ensuring king safety."""
    legal_moves = []
    player_pieces = PLAYER_1_PIECES if player == 1 else PLAYER_2_PIECES
    
    for pos, piece in enumerate(board):
        if piece not in player_pieces:
            continue

        potential_destinations = []
        piece_type = piece[0]

        if piece_type == 'K':
            for d in [-1, 1]:
                if 0 <= pos + d <= 7:
                    potential_destinations.append(pos + d)
        elif piece_type == 'N':
            for d in [-2, 2]:
                if 0 <= pos + d <= 7:
                    potential_destinations.append(pos + d)
        elif piece_type == 'R':
            # Move right
            for i in range(pos + 1, 8):
                potential_destinations.append(i)
                if board[i] != ' ':
                    break
            # Move left
            for i in range(pos - 1, -1, -1):
                potential_destinations.append(i)
                if board[i] != ' ':
                    break
        
        for new_pos in potential_destinations:
            # A move is valid if the destination is not occupied by a friendly piece
            if board[new_pos] not in player_pieces:
                # Simulate the move to check for king safety
                next_board = apply_move(board, (pos, new_pos))
                if not is_in_check(next_board, player):
                    legal_moves.append((pos, new_pos))
    return legal_moves

# --- Game Solving Function (Minimax with Memoization) ---

def solve_game(board, player):
    """
    Recursively solves the game state using minimax.
    Returns: (outcome, ply_to_end)
    - outcome: 'WIN' (for P1), 'LOSE' (for P1), 'DRAW'
    - ply_to_end: number of half-moves (ply) until the outcome
    """
    state_key = (board, player)
    if state_key in memo:
        return memo[state_key]

    # Win condition: Opponent's king is captured
    positions = get_piece_positions(board)
    if 'K2' not in positions: return ('WIN', 0)
    if 'K1' not in positions: return ('LOSE', 0)

    legal_moves = get_legal_moves(board, player)

    # Base case: No legal moves
    if not legal_moves:
        # Checkmate if king is in check, otherwise stalemate
        if is_in_check(board, player):
            result = ('LOSE', 0) if player == 1 else ('WIN', 0)
        else:
            result = ('DRAW', 0)
        memo[state_key] = result
        return result

    child_results = []
    for move in legal_moves:
        next_board = apply_move(board, move)
        outcome, ply = solve_game(next_board, 3 - player)
        child_results.append((outcome, ply))

    # Evaluate results based on current player's goal
    if player == 1:  # Player 1 (Maximizer) wants to WIN, quickly
        wins = [r for r in child_results if r[0] == 'WIN']
        if wins:
            min_ply = min(p for o, p in wins)
            result = ('WIN', 1 + min_ply)
        else:
            draws = [r for r in child_results if r[0] == 'DRAW']
            if draws:
                max_ply = max(p for o, p in draws)
                result = ('DRAW', 1 + max_ply)
            else: # Only losses
                max_ply = max(p for o, p in child_results if o[0] == 'L')
                result = ('LOSE', 1 + max_ply)
    else:  # Player 2 (Minimizer) wants P1 to LOSE, quickly
        losses = [r for r in child_results if r[0] == 'LOSE']
        if losses:
            min_ply = min(p for o, p in losses)
            result = ('LOSE', 1 + min_ply)
        else:
            draws = [r for r in child_results if r[0] == 'DRAW']
            if draws:
                max_ply = max(p for o, p in draws)
                result = ('DRAW', 1 + max_ply)
            else: # Only P1 wins
                max_ply = max(p for o, p in child_results if o[0] == 'W')
                result = ('WIN', 1 + max_ply)

    memo[state_key] = result
    return result

# --- Main Execution ---
def main():
    """
    Calculates the shortest forced win for Player 1.
    """
    initial_board = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
    outcome, ply = solve_game(initial_board, 1)

    if outcome == 'WIN':
        # A win in M ply takes (M+1)//2 turns for the first player.
        # e.g., a win in 3 ply (P1->P2->P1) is on P1's 2nd turn. (3+1)//2 = 2.
        turns = (ply + 1) // 2
        print(turns)
    elif outcome == 'DRAW':
        print("The game is a forced draw. No forced win exists.")
    else: # 'LOSE'
        print("Player 1 cannot force a win.")

if __name__ == '__main__':
    main()