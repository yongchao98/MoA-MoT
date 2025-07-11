import functools

def solve_game():
    """
    Analyzes the game to find the shortest forced win for Player 1.
    """
    memo = {}

    # --- Constants and Board Representation ---
    INITIAL_BOARD = ('K1', 'N1', 'R1', '', '', 'R2', 'N2', 'K2')
    PLAYER_PIECES = {
        1: ('K1', 'N1', 'R1'),
        2: ('K2', 'N2', 'R2')
    }
    PLAYER_KING = {1: 'K1', 2: 'K2'}
    OPPONENT_ROOK = {1: 'R2', 2: 'R1'}
    BOARD_SIZE = 8

    # --- Helper Functions ---
    def get_piece_positions(board_tuple):
        """Gets all piece positions from a board tuple."""
        positions = {}
        for i, piece in enumerate(board_tuple):
            if piece:
                positions[piece] = i
        return positions

    def is_king_in_check(board_tuple, player_num):
        """Checks if a player's king is threatened by the opponent's rook."""
        positions = get_piece_positions(board_tuple)
        king_piece = PLAYER_KING[player_num]
        rook_piece = OPPONENT_ROOK[player_num]

        if king_piece not in positions or rook_piece not in positions:
            return False

        king_pos = positions[king_piece]
        rook_pos = positions[rook_piece]

        start = min(king_pos, rook_pos)
        end = max(king_pos, rook_pos)

        for i in range(start + 1, end):
            if board_tuple[i] != '':
                return False
        return True

    def get_legal_moves(board_tuple, player_num):
        """Generates all legal moves for a given player."""
        legal_moves = []
        positions = get_piece_positions(board_tuple)
        my_pieces = PLAYER_PIECES[player_num]

        for piece, start_pos in positions.items():
            if piece not in my_pieces:
                continue

            piece_type = piece[0]
            potential_ends = []
            if piece_type == 'K':
                if start_pos > 0: potential_ends.append(start_pos - 1)
                if start_pos < BOARD_SIZE - 1: potential_ends.append(start_pos + 1)
            elif piece_type == 'N':
                if start_pos > 1: potential_ends.append(start_pos - 2)
                if start_pos < BOARD_SIZE - 2: potential_ends.append(start_pos + 2)
            elif piece_type == 'R':
                for i in range(start_pos + 1, BOARD_SIZE):
                    potential_ends.append(i)
                    if board_tuple[i] != '': break
                for i in range(start_pos - 1, -1, -1):
                    potential_ends.append(i)
                    if board_tuple[i] != '': break

            for end_pos in potential_ends:
                if board_tuple[end_pos] and board_tuple[end_pos] in my_pieces:
                    continue

                board_list = list(board_tuple)
                board_list[end_pos] = board_list[start_pos]
                board_list[start_pos] = ''
                temp_board = tuple(board_list)

                if not is_king_in_check(temp_board, player_num):
                    legal_moves.append((start_pos, end_pos))

        return legal_moves

    # --- Minimax Solver ---
    def solve(board_tuple, player_num):
        """
        Recursively determines the game outcome from the state (board_tuple, player_num).
        Returns: (outcome, plies) where outcome is 1 for P1 win, 2 for P2 win, 0 for draw.
        """
        state = (board_tuple, player_num)
        if state in memo:
            return memo[state]

        # Base case: Win by king capture
        positions = get_piece_positions(board_tuple)
        if PLAYER_KING[2] not in positions: return (1, 0)
        if PLAYER_KING[1] not in positions: return (2, 0)

        moves = get_legal_moves(board_tuple, player_num)

        # Base case: Checkmate or Stalemate
        if not moves:
            if is_king_in_check(board_tuple, player_num):
                # Checkmate: current player is in check with no moves, opponent wins.
                return (3 - player_num, 0)
            else:
                # Stalemate: no moves, not in check.
                return (0, 0)

        win_plies, loss_plies, draw_plies = [], [], []
        next_player_num = 3 - player_num

        for start_pos, end_pos in moves:
            board_list = list(board_tuple)
            board_list[end_pos] = board_list[start_pos]
            board_list[start_pos] = ''
            outcome, plies = solve(tuple(board_list), next_player_num)

            # Outcomes are always from Player 1's perspective
            if outcome == 1: win_plies.append(plies + 1)
            elif outcome == 2: loss_plies.append(plies + 1)
            else: draw_plies.append(plies + 1)

        # Determine result based on player's goal
        if player_num == 1: # P1 wants to win fast
            if win_plies: result = (1, min(win_plies))
            elif draw_plies: result = (0, min(draw_plies))
            else: result = (2, max(loss_plies)) # Prolong inevitable loss
        else: # P2 wants to win fast
            if loss_plies: result = (2, min(loss_plies)) # P2 win is a P1 loss
            elif draw_plies: result = (0, min(draw_plies))
            else: result = (1, max(win_plies)) # Prolong inevitable loss

        memo[state] = result
        return result

    # --- Main Execution ---
    outcome, plies = solve(INITIAL_BOARD, 1)

    if outcome == 1:
        # Number of turns for Player 1 is ceil(plies / 2)
        turns = (plies + 1) // 2
        print("The first player can force a win in this many turns:")
        print(turns)
    elif outcome == 2:
        print("The second player can force a win.")
    else:
        print("The game ends in a draw if both players play optimally.")

solve_game()