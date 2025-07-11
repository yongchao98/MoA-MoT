import sys
from functools import lru_cache

# Set a higher recursion limit for potentially deep game trees.
sys.setrecursionlimit(5000)

# --- Game Constants ---
P1_PIECES = {'K1', 'N1', 'R1'}
P2_PIECES = {'K2', 'N2', 'R2'}

# --- Helper Functions ---

def is_king_attacked(board, king_player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    king_piece = 'K1' if king_player == 1 else 'K2'
    rook_piece = 'R2' if king_player == 1 else 'R1'

    try:
        king_pos = board.index(king_piece)
    except ValueError:
        return False  # King not on board, can't be attacked

    try:
        rook_pos = board.index(rook_piece)
    except ValueError:
        return False  # Opponent's rook not on board

    start, end = sorted((king_pos, rook_pos))
    is_blocked = any(board[i] != ' ' for i in range(start + 1, end))

    return not is_blocked

def generate_moves(board, player):
    """Generates all legal successor boards for the given player."""
    player_pieces = P1_PIECES if player == 1 else P2_PIECES
    
    for pos, piece in enumerate(board):
        if piece not in player_pieces:
            continue

        # --- King Moves (1 step) ---
        if piece.startswith('K'):
            for move in [-1, 1]:
                new_pos = pos + move
                if 0 <= new_pos < 8 and board[new_pos] not in player_pieces:
                    new_board_list = list(board)
                    new_board_list[pos], new_board_list[new_pos] = ' ', piece
                    new_board = tuple(new_board_list)
                    if not is_king_attacked(new_board, player):
                        yield new_board

        # --- Knight Moves (2 steps) ---
        elif piece.startswith('N'):
            for move in [-2, 2]:
                new_pos = pos + move
                if 0 <= new_pos < 8 and board[new_pos] not in player_pieces:
                    new_board_list = list(board)
                    new_board_list[pos], new_board_list[new_pos] = ' ', piece
                    new_board = tuple(new_board_list)
                    if not is_king_attacked(new_board, player):
                        yield new_board

        # --- Rook Moves (any steps) ---
        elif piece.startswith('R'):
            # Move right
            for new_pos in range(pos + 1, 8):
                if board[new_pos] in player_pieces:
                    break
                new_board_list = list(board)
                new_board_list[pos], new_board_list[new_pos] = ' ', piece
                new_board = tuple(new_board_list)
                if not is_king_attacked(new_board, player):
                    yield new_board
                if board[new_pos] != ' ': # Captured opponent or hit own piece
                    break
            # Move left
            for new_pos in range(pos - 1, -1, -1):
                if board[new_pos] in player_pieces:
                    break
                new_board_list = list(board)
                new_board_list[pos], new_board_list[new_pos] = ' ', piece
                new_board = tuple(new_board_list)
                if not is_king_attacked(new_board, player):
                    yield new_board
                if board[new_pos] != ' ': # Captured opponent or hit own piece
                    break

# --- Minimax Solver with Memoization ---

@lru_cache(maxsize=None)
def solve(board, player):
    """
    Determines the game outcome from the current state using minimax logic.
    Returns: (outcome, num_p1_turns)
      - outcome: 1 for P1 win, -1 for P2 win, 0 for Draw.
      - num_p1_turns: Number of P1 moves to achieve the outcome.
    """
    # Terminal State: Win/Loss by King Capture
    if 'K2' not in board: return (1, 0)
    if 'K1' not in board: return (-1, 0)

    legal_next_boards = list(generate_moves(board, player))

    # Terminal State: Stalemate/Checkmate
    if not legal_next_boards:
        if is_king_attacked(board, player):
            # Checkmate: current player loses
            return (-1, 0) if player == 1 else (1, 0)
        else:
            # Stalemate: draw
            return (0, float('inf'))

    child_results = [solve(next_board, 3 - player) for next_board in legal_next_boards]

    if player == 1:  # Player 1 (maximizer) wants to win fast
        winning_moves = [res for res in child_results if res[0] == 1]
        if winning_moves:
            best_win = min(winning_moves, key=lambda x: x[1])
            return (1, best_win[1] + 1)

        drawing_moves = [res for res in child_results if res[0] == 0]
        if drawing_moves:
            return (0, float('inf'))

        # All moves lead to a loss, so P1 stalls as long as possible
        worst_loss = max(child_results, key=lambda x: x[1])
        return (-1, worst_loss[1])

    else:  # Player 2 (minimizer) wants to win fast or stall a loss
        winning_moves_p2 = [res for res in child_results if res[0] == -1]
        if winning_moves_p2:
            # P2 can win, chooses the fastest path (P1 turns are irrelevant)
            return (-1, 0)

        drawing_moves = [res for res in child_results if res[0] == 0]
        if drawing_moves:
            return (0, float('inf'))

        # All moves lead to a loss for P2, so P2 stalls as long as possible
        best_stall = max(child_results, key=lambda x: x[1])
        return (1, best_stall[1])

def main():
    """
    Calculates and prints the number of turns for Player 1's forced win.
    """
    initial_board = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
    initial_player = 1

    outcome, turns = solve(initial_board, initial_player)

    if outcome == 1:
        print(turns)
    elif outcome == -1:
        print("Player 2 can force a win.")
    else:
        print("The game is a forced draw.")

if __name__ == '__main__':
    main()