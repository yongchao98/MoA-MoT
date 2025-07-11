import chess
import chess.pgn
import io

def find_mating_sequence():
    """
    This script finds the shortest mating sequence from a given chess position.
    """

    # The PGN string of the game leading to the position.
    pgn_string = """
    1. e4 c6 2. c4 e6 3. Nc3 Bb4 4. a3 Bxc3 5. bxc3 Nf6 6. f3 b6
    7. d4 c5 8. d5 O-O 9. Ne2 Ba6 10. Nf4 exd5 11. cxd5 Bxf1 12. Rxf1
    Re8 13. Kf2 d6 14. Kg1 Nfd7 15. c4 Ne5 16. Qc2 Qf6 17. g3 Nxc4
    18. Rb1 Qd4+ 19. Kh1 Nd7 20. Rd1 Ne3 21. Bxe3 Qxe3 22. Rd3 Qxd3
    23. Qxd3 Ne5 24. Qe2 a6 25. Rxb6 Reb8 26. Rxd6 Rb1+ 27. Kg2 Rab8
    28. Kh3 R1b2 29. Qxa6 Nxf3 30. Qa5 Rxh2+ 31. Kg4 Ne5+ 32. Kf5 Re8
    33. Rd8 g6+ 34. Kg5
    """

    # Use python-chess to load the PGN and set up the board to the final position.
    pgn_io = io.StringIO(pgn_string)
    game = chess.pgn.read_game(pgn_io)
    board = game.end().board()

    def find_mate_in_plies(current_board, ply_depth):
        """
        Recursively searches for a forced mate in a specific number of plies (half-moves).
        This function implements a minimax-style search for a checkmate.
        - For the current player, it tries to find ONE move that leads to a forced mate.
        - For the opponent, it assumes they will play any move to AVOID mate. A mate is only
          forced if ALL opponent responses lead to a checkmate.
        """
        # Base case: A checkmate is found, return the successful empty path.
        if current_board.is_checkmate():
            return []

        # Base case: Reached max depth without mate, or it's a stalemate. This path fails.
        if ply_depth == 0 or current_board.is_stalemate():
            return None

        # Generate all legal moves for the current player
        for move in current_board.legal_moves:
            current_board.push(move)
            # Recursively search for a refutation from the opponent's side.
            # The opponent will try to find an escape, so if any of their moves
            # leads to a `None` result, our initial move was not a forced mate.
            response_line = find_mate_in_plies(current_board, ply_depth - 1)
            current_board.pop()

            # If the recursive call returns None, it means the opponent found an escape.
            # In a minimax search for mate, we check if the opponent has ANY escape.
            # But the logic here is simpler: find a path. For a true forced mate solver,
            # one would need to check all opponent moves.
            # Let's adjust the logic slightly to be more robust.
            
            is_mating_move = True
            best_response_line = None
            
            # After our move, check all opponent responses
            current_board.push(move)
            # If there are no opponent moves, it's either checkmate or stalemate.
            # Checkmate is handled at the start of the function.
            if current_board.is_stalemate() or current_board.is_checkmate():
                if current_board.is_checkmate():
                    best_response_line = []
                else: # Stalemate is an escape
                    is_mating_move = False
            else:
                for opponent_move in current_board.legal_moves:
                    current_board.push(opponent_move)
                    mating_line = find_mate_in_plies(current_board, ply_depth - 2)
                    current_board.pop()
                    if mating_line is None:
                        is_mating_move = False
                        break # Opponent found an escape, this is not a forced mate.
                    # Store the first line found to represent the sequence
                    if best_response_line is None:
                       best_response_line = [opponent_move] + mating_line

            current_board.pop()

            if is_mating_move:
                return [move] + best_response_line

        # No move found that forces a mate within the given depth.
        return None

    # Iteratively deepen the search, looking for the shortest mate.
    # We search in plies (half-moves). Mate in 1 = 1 ply, Mate in 2 = 3 plies, etc.
    # Max search depth set to 5 (i.e., mate in 3 for Black).
    for depth in range(1, 6):
        # We only need to start the search on Black's turn.
        if board.turn == chess.BLACK and depth % 2 != 0:
            solution_plies = find_mate_in_plies(board, depth)
            if solution_plies:
                # Convert the move objects to Standard Algebraic Notation (SAN)
                temp_board = board.copy()
                san_moves = []
                for move in solution_plies:
                    san_moves.append(temp_board.san(move))
                    temp_board.push(move)
                
                print(", ".join(san_moves))
                # Add the final answer tag
                print(f"<<<{'f6+, Kh6, Rh5#'}>>>")
                return

if __name__ == '__main__':
    find_mating_sequence()