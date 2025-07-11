import sys
from functools import lru_cache

def solve_game():
    """
    This script solves the described board game to find the minimum number of turns
    for Player 1 to force a win against an optimal opponent.
    """
    # --- Game Setup ---
    BOARD_SIZE = 8
    PLAYER1_PIECES = {"K1", "N1", "R1"}
    PLAYER2_PIECES = {"K2", "N2", "R2"}
    INITIAL_BOARD = ("K1", "N1", "R1", " ", " ", "R2", "N2", "K2")
    INFINITY = float('inf')

    # --- Helper Functions ---
    def get_owner(piece):
        if piece in PLAYER1_PIECES: return 1
        if piece in PLAYER2_PIECES: return 2
        return 0

    def find_piece(board, piece_to_find):
        try:
            return board.index(piece_to_find)
        except ValueError:
            return None

    def apply_move(board, move):
        start_pos, end_pos = move
        piece_to_move = board[start_pos]
        new_board_list = list(board)
        new_board_list[end_pos] = piece_to_move
        new_board_list[start_pos] = " "
        return tuple(new_board_list)

    # --- Core Game Logic ---
    def is_rook_check(board, king_player):
        king_piece = "K1" if king_player == 1 else "K2"
        rook_piece = "R2" if king_player == 1 else "R1"
        king_pos = find_piece(board, king_piece)
        rook_pos = find_piece(board, rook_piece)

        if king_pos is None or rook_pos is None:
            return False

        start, end = min(king_pos, rook_pos), max(king_pos, rook_pos)
        path_clear = all(board[i] == " " for i in range(start + 1, end))
        return path_clear

    def is_move_safe(board, player, move):
        next_board = apply_move(board, move)
        return not is_rook_check(next_board, player)

    def get_legal_moves(board, player):
        legal_moves = []
        player_pieces = PLAYER1_PIECES if player == 1 else PLAYER2_PIECES
        for start_pos, piece in enumerate(board):
            if piece in player_pieces:
                # King Moves (+/- 1)
                if 'K' in piece:
                    for d in [-1, 1]:
                        end_pos = start_pos + d
                        if 0 <= end_pos < BOARD_SIZE and get_owner(board[end_pos]) != player:
                            move = (start_pos, end_pos)
                            if is_move_safe(board, player, move): legal_moves.append(move)
                # Knight Moves (+/- 2)
                elif 'N' in piece:
                    for d in [-2, 2]:
                        end_pos = start_pos + d
                        if 0 <= end_pos < BOARD_SIZE and get_owner(board[end_pos]) != player:
                            move = (start_pos, end_pos)
                            if is_move_safe(board, player, move): legal_moves.append(move)
                # Rook Moves
                elif 'R' in piece:
                    for direction in [-1, 1]:
                        for i in range(1, BOARD_SIZE):
                            end_pos = start_pos + i * direction
                            if not (0 <= end_pos < BOARD_SIZE): break
                            if get_owner(board[end_pos]) == player: break
                            
                            move = (start_pos, end_pos)
                            if is_move_safe(board, player, move): legal_moves.append(move)
                            
                            if get_owner(board[end_pos]) != 0: break # Stop after capture
        return legal_moves

    # --- Minimax Solver ---
    @lru_cache(maxsize=None)
    def find_shortest_win(board, player):
        # Check for terminal states (win/loss by capture)
        if find_piece(board, "K2") is None: return 0
        if find_piece(board, "K1") is None: return INFINITY

        legal_moves = get_legal_moves(board, player)

        # Check for terminal states (checkmate/stalemate)
        if not legal_moves:
            if is_rook_check(board, player): # Checkmate
                return INFINITY if player == 1 else 0
            else: # Stalemate
                return INFINITY

        # Player 1 wants to MINIMIZE the number of moves to win
        if player == 1:
            min_plies = INFINITY
            for move in legal_moves:
                next_board = apply_move(board, move)
                result = find_shortest_win(next_board, 2)
                if result != INFINITY:
                    min_plies = min(min_plies, 1 + result)
            return min_plies
        # Player 2 wants to MAXIMIZE the number of moves to delay losing
        else: # player == 2
            max_plies = -INFINITY
            for move in legal_moves:
                next_board = apply_move(board, move)
                result = find_shortest_win(next_board, 1)
                # If P2 can force a draw/win, they will.
                if result == INFINITY:
                    return INFINITY
                max_plies = max(max_plies, 1 + result)
            
            # If P1 has a forced win, P2 will pick the path that takes the longest.
            return max_plies if max_plies > -INFINITY else INFINITY
            
    # --- Main Execution ---
    # Set a higher recursion limit for deep searches if needed
    # sys.setrecursionlimit(2000)
    
    plies_to_win = find_shortest_win(INITIAL_BOARD, 1)

    if plies_to_win == INFINITY:
        print("Player 1 cannot force a win against an optimal opponent.")
    else:
        # The number of turns for P1 is ceil(plies / 2)
        p1_turns = (plies_to_win + 1) // 2
        
        # Outputting the 'equation' as requested
        print(f"The shortest forced win for Player 1 takes {plies_to_win} plies (half-moves).")
        print(f"The number of turns for Player 1 is calculated as: ({plies_to_win} + 1) // 2 = {p1_turns}")
        print(f"Final Answer: {p1_turns}")
        print(f"<<<{p1_turns}>>>")


solve_game()