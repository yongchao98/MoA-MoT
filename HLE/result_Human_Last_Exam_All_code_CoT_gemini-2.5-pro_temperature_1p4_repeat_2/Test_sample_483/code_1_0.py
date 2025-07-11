import sys

# Setting a higher recursion limit for potentially deep game trees.
sys.setrecursionlimit(2000)

# --- Game Configuration ---
# Piece representations
EMPTY = " "
K1, N1, R1 = "K1", "N1", "R1"
K2, N2, R2 = "K2", "N2", "R2"

PLAYER_1_PIECES = {K1, N1, R1}
PLAYER_2_PIECES = {K2, N2, R2}
INITIAL_BOARD = tuple([K1, N1, R1, EMPTY, EMPTY, R2, N2, K2])
BOARD_SIZE = 8

# --- Solver Configuration ---
# Using integers for outcome types for cleaner code and comparison.
WIN, LOSS, DRAW = 1, -1, 0
# Memoization table to store results of previously computed states.
memo = {}

# --- Helper Functions ---
def get_player(piece):
    """Determines which player a piece belongs to."""
    if piece in PLAYER_1_PIECES:
        return 1
    if piece in PLAYER_2_PIECES:
        return 2
    return None

def is_king_in_check(board, player):
    """Checks if the specified player's King is under attack by an opponent's Rook."""
    king_piece = K1 if player == 1 else K2
    rook_piece = R2 if player == 1 else R1

    try:
        king_pos = board.index(king_piece)
    except ValueError:
        # If the king is not on the board, it has been captured.
        # This is a terminal state that results in a loss.
        return True

    try:
        rook_pos = board.index(rook_piece)
    except ValueError:
        # If the opponent's rook is not on the board, there is no threat from it.
        return False

    # Check for any blocking pieces between the King and the Rook.
    start, end = sorted((king_pos, rook_pos))
    for i in range(start + 1, end):
        if board[i] != EMPTY:
            return False # The path is blocked.
    # If the path is clear, the King is in check.
    return True

# --- Move Generation Functions ---
def apply_move(board, move):
    """Applies a move to the board and returns the new board state as a tuple."""
    from_pos, to_pos = move
    new_board = list(board)
    piece = new_board[from_pos]
    new_board[to_pos] = piece
    new_board[from_pos] = EMPTY
    return tuple(new_board)

def generate_legal_moves(board, player):
    """Generates all legal moves for a given player, filtering for King safety."""
    player_pieces = PLAYER_1_PIECES if player == 1 else PLAYER_2_PIECES
    
    potential_moves = []
    # 1. Generate all possible moves based on piece movement rules.
    for pos, piece in enumerate(board):
        if piece not in player_pieces:
            continue

        if piece in {K1, K2}: # King movement
            for d in [-1, 1]:
                new_pos = pos + d
                if 0 <= new_pos < BOARD_SIZE and get_player(board[new_pos]) != player:
                    potential_moves.append((pos, new_pos))
        elif piece in {N1, N2}: # Knight movement
            for d in [-2, 2]:
                new_pos = pos + d
                if 0 <= new_pos < BOARD_SIZE and get_player(board[new_pos]) != player:
                    potential_moves.append((pos, new_pos))
        elif piece in {R1, R2}: # Rook movement
            # Move left
            for new_pos in range(pos - 1, -1, -1):
                if get_player(board[new_pos]) == player: break
                potential_moves.append((pos, new_pos))
                if get_player(board[new_pos]) is not None: break
            # Move right
            for new_pos in range(pos + 1, BOARD_SIZE):
                if get_player(board[new_pos]) == player: break
                potential_moves.append((pos, new_pos))
                if get_player(board[new_pos]) is not None: break

    # 2. Filter out moves that leave the player's own King in check.
    legal_moves = []
    for move in potential_moves:
        temp_board = apply_move(board, move)
        if not is_king_in_check(temp_board, player):
            legal_moves.append(move)
            
    return legal_moves

# --- Minimax Solver ---
def solve(board, player):
    """
    Recursively solves the game from the given state using minimax with memoization.
    Returns a tuple: (outcome, ply), where outcome is WIN/LOSS/DRAW and ply
    is the number of half-moves to that outcome.
    """
    state = (board, player)
    if state in memo:
        return memo[state]

    opponent = 2 if player == 1 else 1
    opponent_king = K2 if player == 1 else K1

    legal_moves = generate_legal_moves(board, player)

    # Base Case 1: Terminal state (no legal moves).
    if not legal_moves:
        if is_king_in_check(board, player):
             # Checkmated: The current player loses immediately.
            result = (LOSS, 0)
        else:
             # Stalemate: The game is a draw.
            result = (DRAW, 0)
        memo[state] = result
        return result

    possible_outcomes = []
    for move in legal_moves:
        new_board = apply_move(board, move)

        # Base Case 2: Immediate win by capturing the opponent's king.
        if opponent_king not in new_board:
            possible_outcomes.append((WIN, 1))
            continue
        
        # Recursive Step: Solve for the opponent's turn.
        res_type, res_ply = solve(new_board, opponent)

        # Translate the opponent's result to the current player's perspective.
        if res_type == WIN:
            outcome = (LOSS, res_ply + 1) # If opponent wins, I lose.
        elif res_type == LOSS:
            outcome = (WIN, res_ply + 1) # If opponent loses, I win.
        else: # DRAW
            outcome = (DRAW, res_ply + 1)
        possible_outcomes.append(outcome)

    # Aggregate results: Choose the best move according to minimax principles.
    wins = [o for o in possible_outcomes if o[0] == WIN]
    if wins:
        # If there are winning moves, choose the one with the fewest ply (fastest win).
        min_ply = min(w[1] for w in wins)
        best_outcome = (WIN, min_ply)
    else:
        draws = [o for o in possible_outcomes if o[0] == DRAW]
        if draws:
            # If no winning moves, but can force a draw, do so.
            best_outcome = (DRAW, 0)
        else:
            # If all moves lead to a loss, choose the move that stalls the longest (slowest loss).
            losses = [o for o in possible_outcomes if o[0] == LOSS]
            max_ply = max(l[1] for l in losses)
            best_outcome = (LOSS, max_ply)
            
    memo[state] = best_outcome
    return best_outcome

# --- Main Execution Block ---
def main():
    """
    Runs the solver and prints the result.
    """
    # Start the search from the initial state for Player 1.
    outcome, ply = solve(INITIAL_BOARD, 1)

    if outcome == WIN:
        # A "turn" consists of a move by Player 1 and a move by Player 2.
        # "ply" is the number of half-moves. Player 1 wins on an odd-numbered ply.
        # A win on ply 'd' happens on turn number '(d+1)/2'.
        # For example, a win in 1 ply is on turn 1. A win in 3 ply is on turn 2.
        turns = (ply + 1) // 2
        print(f"The first player can force a win.")
        print(f"The minimum number of ply (half-moves) to guarantee a win is: {ply}")
        print(f"The final number of turns is derived from the number of ply.")
        print(f"Calculation: Turns = (Ply + 1) // 2")
        print(f"Equation: {turns} = ({ply} + 1) // 2")
        print(f"Final Answer: It would take {turns} turns for the first player to force a win.")
    elif outcome == DRAW:
        print("The game ends in a forced draw. Player 1 cannot force a win.")
    else: # LOSS
        print("Player 1 is in a losing position and cannot force a win.")

if __name__ == "__main__":
    main()
