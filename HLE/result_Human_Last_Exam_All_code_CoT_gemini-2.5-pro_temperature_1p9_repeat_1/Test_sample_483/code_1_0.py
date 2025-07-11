import sys

# Setting a higher recursion limit for the recursive solver, a good practice.
sys.setrecursionlimit(2000)

# --- Game Configuration ---
INITIAL_BOARD = ('K1', 'N1', 'R1', ' ', ' ', 'R2', 'N2', 'K2')
P1_PIECES = {'K1', 'N1', 'R1'}
P2_PIECES = {'K2', 'N2', 'R2'}

# --- Solver Constants ---
WIN = 1
DRAW = 0
LOSS = -1

# Memoization cache to store results of solved states
memo = {}

# --- Helper Functions ---

def get_player_info(player):
    """Returns the pieces belonging to the current player and the opponent."""
    if player == 1:
        return P1_PIECES, P2_PIECES
    else:
        return P2_PIECES, P1_PIECES

def is_king_attacked(board, player):
    """Checks if the specified player's king is under attack by the opponent's rook."""
    try:
        king_piece = 'K1' if player == 1 else 'K2'
        rook_piece = 'R2' if player == 1 else 'R1'
        
        king_pos = board.index(king_piece)
        rook_pos = board.index(rook_piece)
    except ValueError:
        # A key piece is not on the board, so no attack is possible.
        return False

    start = min(king_pos, rook_pos) + 1
    end = max(king_pos, rook_pos)
    
    # Check for blocking pieces between the king and rook.
    for i in range(start, end):
        if board[i] != ' ':
            return False  # Path is blocked.
    
    return True  # Path is clear, king is attacked.

def get_legal_moves(board, player):
    """Generates all legal moves for a given player."""
    legal_moves = []
    player_pieces, opponent_pieces = get_player_info(player)

    for start_pos, piece in enumerate(board):
        if piece in player_pieces:
            piece_type = piece[0]
            
            potential_dests = []
            if piece_type == 'K':  # King moves
                potential_dests.extend([start_pos - 1, start_pos + 1])
            elif piece_type == 'N':  # Knight moves
                potential_dests.extend([start_pos - 2, start_pos + 2])
            elif piece_type == 'R':  # Rook moves
                # Scan left
                for i in range(start_pos - 1, -1, -1):
                    if board[i] == ' ':
                        potential_dests.append(i)
                    else:
                        if board[i] in opponent_pieces:
                            potential_dests.append(i)
                        break
                # Scan right
                for i in range(start_pos + 1, 8):
                    if board[i] == ' ':
                        potential_dests.append(i)
                    else:
                        if board[i] in opponent_pieces:
                            potential_dests.append(i)
                        break
            
            for dest_pos in potential_dests:
                # Validate destination: bounds and not occupied by a friendly piece
                if not (0 <= dest_pos < 8) or board[dest_pos] in player_pieces:
                    continue
                
                # Simulate the move and check for self-check
                board_list = list(board)
                board_list[dest_pos] = piece
                board_list[start_pos] = ' '
                temp_board = tuple(board_list)
                
                if not is_king_attacked(temp_board, player):
                    legal_moves.append((start_pos, dest_pos))
                    
    return legal_moves

def apply_move(board, move):
    """Applies a move and returns the new board state and the captured piece."""
    start_pos, end_pos = move
    board_list = list(board)
    
    piece = board_list[start_pos]
    captured_piece = board_list[end_pos]
    
    board_list[end_pos] = piece
    board_list[start_pos] = ' '
    
    return tuple(board_list), captured_piece

def solve(board, player):
    """
    Recursively solves the game state using minimax with memoization.
    Returns a tuple (outcome, moves) for the current player.
    """
    state = (board, player)
    if state in memo:
        return memo[state]

    opponent = 2 if player == 1 else 1
    opponent_king = 'K2' if player == 1 else 'K1'
    
    legal_moves = get_legal_moves(board, player)

    # Base case: No legal moves available
    if not legal_moves:
        if is_king_attacked(board, player):
            result = (LOSS, 0)  # Checkmated
        else:
            result = (DRAW, 0)  # Stalemate
        memo[state] = result
        return result

    outcomes = []
    for move in legal_moves:
        new_board, captured_piece = apply_move(board, move)
        
        # Check for an immediate win by capturing the opponent's king
        if captured_piece == opponent_king:
            outcomes.append((WIN, 1))
        else:
            # Recurse for the opponent's turn
            opp_outcome, opp_moves = solve(new_board, opponent)
            
            # Translate opponent's outcome to the current player's perspective
            if opp_outcome == WIN:
                my_outcome = LOSS
            elif opp_outcome == LOSS:
                my_outcome = WIN
            else: # DRAW
                my_outcome = DRAW
            
            my_moves = opp_moves + 1 if my_outcome != DRAW else 0
            outcomes.append((my_outcome, my_moves))

    # Determine the best outcome based on minimax principles
    wins = [o for o in outcomes if o[0] == WIN]
    draws = [o for o in outcomes if o[0] == DRAW]
    
    if wins:
        # If winning is possible, choose the fastest win
        min_moves = min(m for _, m in wins)
        result = (WIN, min_moves)
    elif draws:
        result = (DRAW, 0)
    else: # Only losses are available
        # If losing is inevitable, choose the slowest loss (stall)
        losses = [o for o in outcomes if o[0] == LOSS]
        max_moves = max(m for _, m in losses)
        result = (LOSS, max_moves)

    memo[state] = result
    return result

if __name__ == "__main__":
    # Start the analysis from the initial board state for Player 1
    outcome, moves = solve(INITIAL_BOARD, 1)

    if outcome == WIN:
        print(moves)
    elif outcome == DRAW:
        print("The game ends in a forced draw.")
    else: # LOSS
        print(f"Player 1 has a forced loss in {moves} moves.")
