import sys

def solve():
    """
    This script solves the described board game problem by finding the
    shortest forced win for Player 1 using a minimax search algorithm.
    """

    # It's a search-heavy problem, so increasing the recursion limit is a safe measure.
    sys.setrecursionlimit(2000)

    # Memoization table to cache results of previously seen game states
    MEMO = {}

    # Piece sets for easy identification of a piece's owner
    P1_PIECES = {'K1', 'N1', 'R1'}
    P2_PIECES = {'K2', 'N2', 'R2'}

    def get_piece_positions(board):
        """Helper function to get a map of piece names to their board positions."""
        positions = {}
        for i, piece in enumerate(board):
            if piece:
                positions[piece] = i
        return positions

    def is_king_in_check(board, player):
        """
        Checks if the specified player's king is under attack by an opponent's rook.
        A king is under attack if an opponent's rook has a clear line of sight
        with no other pieces in between.
        """
        positions = get_piece_positions(board)
        king_piece = 'K1' if player == 1 else 'K2'
        opponent_rook_piece = 'R2' if player == 1 else 'R1'

        # If king or opponent's rook is not on the board, no check is possible.
        if king_piece not in positions or opponent_rook_piece not in positions:
            return False

        king_pos = positions[king_piece]
        rook_pos = positions[opponent_rook_piece]

        start = min(king_pos, rook_pos)
        end = max(king_pos, rook_pos)

        # Check for any piece on the squares between the king and the rook.
        for i in range(start + 1, end):
            if board[i] != '':
                return False  # Path is blocked, king is safe.

        return True  # Path is clear, king is in check.

    def get_legal_moves(board, player):
        """
        Generates all legal next board states for the given player.
        """
        legal_next_boards = set()
        player_pieces = P1_PIECES if player == 1 else P2_PIECES

        for i, piece in enumerate(board):
            if piece in player_pieces:
                piece_type = piece[0]
                potential_destinations = []

                if piece_type == 'K':  # King moves
                    potential_destinations.extend([i - 1, i + 1])
                elif piece_type == 'N':  # Knight moves
                    potential_destinations.extend([i - 2, i + 2])
                elif piece_type == 'R':  # Rook moves
                    # Move right
                    for j in range(i + 1, 8):
                        potential_destinations.append(j)
                        if board[j] != '':
                            break
                    # Move left
                    for j in range(i - 1, -1, -1):
                        potential_destinations.append(j)
                        if board[j] != '':
                            break
                
                for dest in potential_destinations:
                    # Check if move is within board boundaries and not on a friendly piece
                    if 0 <= dest < 8 and board[dest] not in player_pieces:
                        # Simulate the move
                        new_board_list = list(board)
                        new_board_list[dest] = piece
                        new_board_list[i] = ''
                        new_board = tuple(new_board_list)

                        # A move is only legal if it doesn't leave the king in check
                        if not is_king_in_check(new_board, player):
                            legal_next_boards.add(new_board)

        return list(legal_next_boards)

    def solve_game_state(board, player):
        """
        Recursively determines the outcome of a game state using the minimax algorithm.
        Returns a tuple (outcome, plies_to_outcome).
        """
        state = (board, player)
        if state in MEMO:
            return MEMO[state]

        # Base case: A king has been captured (win/loss)
        if 'K2' not in board:
            return (1, 0)  # P1 wins
        if 'K1' not in board:
            return (2, 0)  # P2 wins

        # Get all legal moves for the current player
        legal_moves = get_legal_moves(board, player)

        # Base case: No legal moves available (checkmate or stalemate)
        if not legal_moves:
            if is_king_in_check(board, player):
                # Checkmate: Current player loses
                result = (2 if player == 1 else 1, 0)
            else:
                # Stalemate: Draw
                result = (0, 0)
            MEMO[state] = result
            return result

        # Recursive step: Evaluate outcomes of all possible next states
        next_player = 2 if player == 1 else 1
        child_outcomes = [solve_game_state(move, next_player) for move in legal_moves]

        if player == 1:  # P1 wants to win as quickly as possible
            p1_wins = [res for res in child_outcomes if res[0] == 1]
            if p1_wins:
                min_plies = min(res[1] for res in p1_wins)
                result = (1, min_plies + 1)
            else:
                draws = [res for res in child_outcomes if res[0] == 0]
                if draws:
                    max_plies = max(res[1] for res in draws)
                    result = (0, max_plies + 1)
                else: # Must lose, so stall as long as possible
                    losses = [res for res in child_outcomes if res[0] == 2]
                    max_plies = max(res[1] for res in losses)
                    result = (2, max_plies + 1)
        
        else:  # P2 wants to win or stall as long as possible
            p2_wins = [res for res in child_outcomes if res[0] == 2]
            if p2_wins:
                min_plies = min(res[1] for res in p2_wins)
                result = (2, min_plies + 1)
            else:
                draws = [res for res in child_outcomes if res[0] == 0]
                if draws:
                    max_plies = max(res[1] for res in draws)
                    result = (0, max_plies + 1)
                else: # Must lose, so stall as long as possible
                    losses = [res for res in child_outcomes if res[0] == 1]
                    max_plies = max(res[1] for res in losses)
                    result = (1, max_plies + 1)

        MEMO[state] = result
        return result

    # --- Main Execution ---
    initial_board = ('K1', 'N1', 'R1', '', '', 'R2', 'N2', 'K2')
    final_outcome, total_plies = solve_game_state(initial_board, 1)

    # The number of turns for P1 is the ceiling of half the total plies.
    if final_outcome == 1:
        p1_turns = (total_plies + 1) // 2
        print(p1_turns)
    else:
        # This case implies P1 cannot force a win.
        # Based on the problem phrasing, a forced win is expected.
        print("Player 1 cannot force a win.")

solve()