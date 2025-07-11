import chess

def solve_chess_checkmates():
    """
    Calculates the number of legal chess positions where White (King, Pawn, Knight)
    has checkmated Black (King).

    This is a computationally intensive task and may take several minutes to run.
    """
    count = 0
    all_squares = list(range(64))
    
    # To be a valid position, a white pawn cannot be on the 1st or 8th rank.
    # This filters the squares a pawn can occupy.
    pawn_squares = [sq for sq in all_squares if chess.square_rank(sq) not in [0, 7]]

    print("Calculating all possible checkmate positions. This may take several minutes...")

    # Iterate through every possible square for the Black King (k).
    for k_sq in all_squares:
        # Iterate through every possible square for the White King (K).
        for K_sq in all_squares:
            if K_sq == k_sq:
                continue

            # Legality check 1: The kings cannot be on adjacent squares.
            if chess.square_distance(k_sq, K_sq) <= 1:
                continue

            # Iterate through every possible square for the White Knight (N).
            for N_sq in all_squares:
                if N_sq == k_sq or N_sq == K_sq:
                    continue
                
                # Iterate through the valid squares for the White Pawn (P).
                for p_sq in pawn_squares:
                    if p_sq == k_sq or p_sq == K_sq or p_sq == N_sq:
                        continue

                    # At this point, we have a unique placement for all 4 pieces.
                    # Now we set up this position on a board to check if it's a mate.
                    board = chess.Board(fen=None)  # Create an empty board.
                    board.set_piece_at(k_sq, chess.Piece(chess.KING, chess.BLACK))
                    board.set_piece_at(K_sq, chess.Piece(chess.KING, chess.WHITE))
                    board.set_piece_at(N_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
                    board.set_piece_at(p_sq, chess.Piece(chess.PAWN, chess.WHITE))

                    # For the position to be a checkmate, it must be the turn
                    # of the player who is mated.
                    board.turn = chess.BLACK

                    # The is_checkmate() method efficiently checks if the current player
                    # (Black) is in check and has no legal moves.
                    # Our prior checks ensure the position's basic legality.
                    if board.is_checkmate():
                        count += 1
                        
    print("The total number of checkmate positions is:")
    print(count)

if __name__ == '__main__':
    solve_chess_checkmates()