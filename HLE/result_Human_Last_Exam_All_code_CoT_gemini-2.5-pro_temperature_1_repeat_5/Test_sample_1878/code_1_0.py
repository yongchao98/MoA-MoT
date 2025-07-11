# To run this code, you first need to install the `python-chess` library.
# You can do this by running the following command in your terminal or command prompt:
# pip install python-chess

import chess
import itertools

def find_checkmate_positions():
    """
    This function systematically checks all possible board positions
    to find how many are legal checkmates for White, who has a King,
    a Knight, and a Pawn, against Black's lone King.
    """
    mate_count = 0
    squares = range(64)

    # We use itertools.permutations to generate every possible arrangement of
    # the 4 pieces on the 64 squares of the chessboard.
    # The pieces are: White King, White Knight, White Pawn, Black King.
    for p in itertools.permutations(squares, 4):
        wk_sq, wn_sq, wp_sq, bk_sq = p

        # Condition 1: The White Pawn cannot be on its starting or promotion rank.
        # A pawn on these ranks is illegal (it should have promoted on rank 8).
        # chess.square_rank() returns 0 for rank '1' and 7 for rank '8'.
        wp_rank = chess.square_rank(wp_sq)
        if wp_rank == 0 or wp_rank == 7:
            continue

        # Condition 2: The two Kings cannot be on adjacent squares.
        # This is an illegal position in chess. chess.distance() calculates
        # the king-move distance between two squares.
        if chess.distance(wk_sq, bk_sq) <= 1:
            continue

        # If the basic placement is valid, create a board for deeper analysis.
        board = chess.Board(fen=None)  # Create a blank board
        board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(wn_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_sq, chess.Piece(chess.PAWN, chess.WHITE))
        board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))

        # Condition 3: The position must be "legal". A simple and crucial check
        # is that the player whose turn it isn't cannot be in check.
        # Here, we check if the White King is attacked by the Black King.
        board.turn = chess.WHITE
        if board.is_check():
            continue

        # Condition 4: The position must be checkmate for Black.
        # We set the turn to Black and use the library's is_checkmate() method.
        board.turn = chess.BLACK
        if board.is_checkmate():
            mate_count += 1
    
    # The final count is the result. This print statement outputs the number.
    print(f"The total number of constructible checkmate positions is: {mate_count}")

# Execute the function to find and print the result.
# Note: This script will take a few minutes to run due to the large number of positions to check.
find_checkmate_positions()