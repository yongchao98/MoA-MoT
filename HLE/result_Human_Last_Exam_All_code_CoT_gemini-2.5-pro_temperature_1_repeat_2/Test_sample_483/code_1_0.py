import sys

# Increase recursion limit for deep game trees.
sys.setrecursionlimit(2000)

# Constants
BOARD_SIZE = 8
P1_PIECES = {'K1', 'N1', 'R1'}
P2_PIECES = {'K2', 'N2', 'R2'}
P1_KING = 'K1'
P2_KING = 'K2'
P1_ROOK = 'R1'
P2_ROOK = 'R2'

# Memoization table to store results of solved states
memo = {}

# Outcome constants
WIN = 1
LOSS = -1
DRAW = 0

def is_king_in_check(board, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    king_piece = P1_KING if player == 'P1' else P2_KING
    opponent_rook_piece = P2_ROOK if player == 'P1' else P1_ROOK

    try:
        king_pos = board.index(king_piece)
    except ValueError:
        return False  # King is captured, so not in check

    try:
        rook_pos = board.index(opponent_rook_piece)
    except ValueError:
        return False  # Opponent's rook is captured, no threat

    start = min(king_pos, rook_pos) + 1
    end = max(king_pos, rook_pos)

    for i in range(start, end):
        if board[i] is not None:
            return False  # Path is blocked, king is safe

    return True  # Path is clear, king is in check

def generate_moves(board, player):
    """Generates all legal moves for a given player."""
    legal_moves = []
    player_pieces = P1_PIECES if player == 'P1' else P2_PIECES
    opponent_pieces = P2_PIECES if player == 'P1' else P1_PIECES

    for pos, piece in enumerate(board):
        if piece in player_pieces:
            piece_type = piece[0]

            dests = []
            if piece_type == 'K':  # King
                dests = [pos - 1, pos + 1]
            elif piece_type == 'N':  # Knight
                dests = [pos - 2, pos + 2]
            elif piece_type == 'R':  # Rook
                # Move right
                for i in range(pos + 1, BOARD_SIZE):
                    dests.append(i)
                    if board[i] is not None:
                        break
                # Move left
                for i in range(pos - 1, -1, -1):
                    dests.append(i)
                    if board[i] is not None:
                        break

            for dest in dests:
                if 0 <= dest < BOARD_SIZE and (board[dest] is None or board[dest] in opponent_pieces):
                    new_board_list = list(board)
                    new_board_list[dest] = piece
                    new_board_list[pos] = None
                    new_board = tuple(new_board_list)

                    if not is_king_in_check(new_board, player):
                        legal_moves.append(new_board)
    return legal_moves

def solve(board, player):
    """
    Recursively solves the game using minimax.
    Returns a tuple (outcome, plies).
    """
    state_key = (board, player)
    if state_key in memo:
        return memo[state_key]

    if P2_KING not in board: return (WIN, 0)
    if P1_KING not in board: return (LOSS, 0)

    legal_next_boards = generate_moves(board, player)

    if not legal_next_boards:
        if is_king_in_check(board, player):
            return (LOSS, 0) if player == 'P1' else (WIN, 0)
        else:
            return (DRAW, 0)

    outcomes = []
    next_player = 'P2' if player == 'P1' else 'P1'
    for next_board in legal_next_boards:
        outcome, plies = solve(next_board, next_player)
        outcomes.append((outcome, plies + 1))

    if player == 'P1':  # Maximizer (P1 wants to win fast)
        wins = [r for r in outcomes if r[0] == WIN]
        if wins:
            result = min(wins, key=lambda x: x[1])
        else:
            draws = [r for r in outcomes if r[0] == DRAW]
            if draws:
                result = max(draws, key=lambda x: x[1])
            else:
                result = max(outcomes, key=lambda x: x[1])
    else:  # Minimizer (P2 wants to lose slow)
        losses = [r for r in outcomes if r[0] == LOSS]
        if losses:
            result = min(losses, key=lambda x: x[1])
        else:
            draws = [r for r in outcomes if r[0] == DRAW]
            if draws:
                result = max(draws, key=lambda x: x[1])
            else:
                result = max(outcomes, key=lambda x: x[1])

    memo[state_key] = result
    return result

def main():
    """
    Sets up the initial board and starts the solver.
    """
    initial_board = (
        'K1', 'N1', 'R1', None, None, 'R2', 'N2', 'K2'
    )
    
    # Start the search from Player 1's turn.
    outcome, plies = solve(initial_board, 'P1')

    if outcome == WIN:
        # A "turn" consists of one move by P1 and potentially one by P2.
        # plies=1 (P1 wins on first move) -> Turn 1
        # plies=3 (P1, P2, P1 wins) -> Turn 2
        # So, turns = ceil(plies / 2)
        turns = (plies + 1) // 2
        print(turns)
        print(f"<<<{turns}>>>")
    elif outcome == DRAW:
        print("DRAW")
        print("<<<DRAW>>>")
    else: # LOSS
        print("LOSS")
        print("<<<LOSS>>>")

if __name__ == '__main__':
    main()