import sys

def solve_game():
    """
    This script solves the described 1D-chess game to find the minimum
    number of turns for Player 1 to force a checkmate.
    """

    # Set a higher recursion limit for potentially deep game tree searches.
    # The default limit of 1000 should be sufficient, but this is a safeguard.
    sys.setrecursionlimit(2000)

    # memo is a cache to store results for states we've already solved.
    memo = {}
    WIN_P1 = "P1_WINS"
    WIN_P2 = "P2_WINS"
    DRAW = "DRAW"
    INFINITY = float('inf')

    def get_piece_player(piece_name: str) -> int | None:
        """Gets the player number (1 or 2) from a piece string like 'K1'."""
        if not piece_name or piece_name.strip() == '':
            return None
        return int(piece_name[1])

    def find_piece(board: tuple, piece_name: str) -> int:
        """Finds the position of a piece on the board."""
        try:
            return board.index(piece_name)
        except ValueError:
            return -1  # Piece not on board (captured)

    def is_king_in_check(board: tuple, player: int) -> bool:
        """Checks if the specified player's king is under attack by the opponent's rook."""
        king_piece = f'K{player}'
        opponent_rook = f'R{3 - player}'

        king_pos = find_piece(board, king_piece)
        rook_pos = find_piece(board, opponent_rook)

        if king_pos == -1 or rook_pos == -1:
            return False

        # Check for a clear line of sight between the king and the rook.
        start, end = sorted((king_pos, rook_pos))
        for i in range(start + 1, end):
            if board[i].strip() != '':
                return False  # Path is blocked by another piece.
        return True

    def generate_legal_moves(board: tuple, player: int) -> list[tuple]:
        """Generates all legal next board states for a given player."""
        legal_next_boards = []
        player_piece_positions = [i for i, p in enumerate(board) if p.strip() and get_piece_player(p) == player]

        for pos in player_piece_positions:
            piece = board[pos]
            piece_type = piece[0]
            
            candidate_dests = []
            # Generate candidate destination squares based on piece type.
            if piece_type == 'K':
                candidate_dests.extend([pos - 1, pos + 1])
            elif piece_type == 'N':
                candidate_dests.extend([pos - 2, pos + 2])
            elif piece_type == 'R':
                # Move right
                for i in range(pos + 1, 8):
                    if board[i].strip() == '':
                        candidate_dests.append(i)
                    else: # A piece is encountered
                        if get_piece_player(board[i]) != player:
                            candidate_dests.append(i) # Can capture opponent
                        break # Rook is blocked
                # Move left
                for i in range(pos - 1, -1, -1):
                    if board[i].strip() == '':
                        candidate_dests.append(i)
                    else: # A piece is encountered
                        if get_piece_player(board[i]) != player:
                            candidate_dests.append(i) # Can capture opponent
                        break # Rook is blocked
            
            # Filter candidates to find strictly legal moves.
            for dest in candidate_dests:
                # Rule 1: Move must be on the board.
                if not (0 <= dest <= 7):
                    continue
                
                # Rule 2: Cannot capture your own piece (already handled by Rook logic).
                dest_piece = board[dest]
                if dest_piece.strip() != '' and get_piece_player(dest_piece) == player:
                    continue
                
                # Create the potential next board state.
                temp_board_list = list(board)
                temp_board_list[dest] = piece
                temp_board_list[pos] = '  '
                next_board = tuple(temp_board_list)
                
                # Rule 3: The move is only legal if the player's king is NOT in check after the move.
                if not is_king_in_check(next_board, player):
                    legal_next_boards.append(next_board)
                    
        return legal_next_boards

    def solve(state: tuple) -> tuple:
        """
        Recursively solves the game from a given state using minimax with memoization.
        A state is a tuple: (board_tuple, player_to_move)
        Returns a tuple: (outcome_string, plies_to_outcome)
        """
        if state in memo:
            return memo[state]

        board, player = state
        opponent = 3 - player
        
        legal_next_boards = generate_legal_moves(board, player)

        # Base Case: Terminal State (No Legal Moves).
        if not legal_next_boards:
            if is_king_in_check(board, player):
                # Checkmate: The current player has lost.
                result = (f'P{opponent}_WINS', 0)
            else:
                # Stalemate: The game is a draw.
                result = (DRAW, 0)
            memo[state] = result
            return result

        # Recursive Step: Explore child states.
        child_results = [solve((next_board, opponent)) for next_board in legal_next_boards]
        
        # Minimax Logic: Determine best outcome based on player's goals.
        if player == 1:  # Player 1 wants to win as fast as possible.
            p1_wins = [r for r in child_results if r[0] == WIN_P1]
            if p1_wins:
                min_plies = min(r[1] for r in p1_wins)
                result = (WIN_P1, 1 + min_plies)
            else:
                draws = [r for r in child_results if r[0] == DRAW]
                if draws:
                    result = (DRAW, INFINITY) # Choose draw over a loss.
                else: # All moves lead to a loss.
                    max_plies = max(r[1] for r in child_results) # Delay loss.
                    result = (WIN_P2, 1 + max_plies)
        
        else:  # Player 2 wants to win, or draw, or prolong loss.
            p2_wins = [r for r in child_results if r[0] == WIN_P2]
            if p2_wins:
                min_plies = min(r[1] for r in p2_wins)
                result = (WIN_P2, 1 + min_plies)
            else:
                draws = [r for r in child_results if r[0] == DRAW]
                if draws:
                    result = (DRAW, INFINITY) # Choose draw over a loss.
                else: # Loss to P1 is inevitable.
                    max_plies = max(r[1] for r in child_results) # Stall for as long as possible.
                    result = (WIN_P1, 1 + max_plies)

        memo[state] = result
        return result

    # --- Main Execution ---
    initial_board = ('K1', 'N1', 'R1', '  ', '  ', 'R2', 'N2', 'K2')
    initial_state = (initial_board, 1)  # Player 1 starts.

    outcome, plies = solve(initial_state)

    if outcome == WIN_P1:
        # A ply is a single move. A turn involves a move by each player.
        # If P1 forces a win in 'k' plies, P1 makes ceil(k/2) moves.
        # This can be calculated with integer division: (plies + 1) // 2.
        turns = (plies + 1) // 2
        print(turns)
    elif outcome == WIN_P2:
        print("Player 1 cannot force a win. Player 2 has a winning strategy.")
    else:
        print("Player 1 cannot force a win. The game ends in a draw with optimal play.")

solve_game()