# This script requires the python-chess library.
# You can install it by running: pip install python-chess

import chess

def count_checkmate_positions():
    """
    Calculates the number of legal checkmate positions with White (K, N, P) vs Black (k).
    
    A position is considered "legal" and a "checkmate" if:
    1. All four pieces are on unique squares.
    2. The White Pawn is on ranks 2-7.
    3. The White and Black kings are not on adjacent squares.
    4. It is Black's turn to move, and Black is in checkmate.
    """
    checkmate_count = 0
    all_squares = chess.SQUARES
    # Pawns cannot be on the 1st or 8th rank (squares with rank index 0 or 7)
    pawn_squares = [sq for sq in all_squares if chess.square_rank(sq) not in [0, 7]]

    # This is a large computation, iterating through millions of positions.
    # It may take a few minutes to run.

    # Iterate over all possible squares for the White King
    for wk_square in all_squares:
        # Iterate over all possible squares for the White Knight
        for wn_square in all_squares:
            if wn_square == wk_square:
                continue

            # Iterate over valid squares for the White Pawn
            for wp_square in pawn_squares:
                if wp_square in [wk_square, wn_square]:
                    continue

                # Iterate over all possible squares for the Black King
                for bk_square in all_squares:
                    if bk_square in [wk_square, wn_square, wp_square]:
                        continue
                        
                    # Basic legality check: Kings cannot be adjacent.
                    # The distance between king squares must be greater than 1.
                    if chess.square_distance(wk_square, bk_square) <= 1:
                        continue

                    # Set up the board position
                    board = chess.Board(fen=None)  # Create an empty board
                    board.set_piece_at(wk_square, chess.Piece(chess.KING, chess.WHITE))
                    board.set_piece_at(wn_square, chess.Piece(chess.KNIGHT, chess.WHITE))
                    board.set_piece_at(wp_square, chess.Piece(chess.PAWN, chess.WHITE))
                    board.set_piece_at(bk_square, chess.Piece(chess.KING, chess.BLACK))

                    # Set the turn to Black to check if Black is mated
                    board.turn = chess.BLACK

                    # The is_checkmate() function checks if the current side to move
                    # is in check and has no legal moves.
                    if board.is_checkmate():
                        checkmate_count += 1
                        
    print(f"The total number of constructible checkmate positions is: {checkmate_count}")

if __name__ == '__main__':
    count_checkmate_positions()