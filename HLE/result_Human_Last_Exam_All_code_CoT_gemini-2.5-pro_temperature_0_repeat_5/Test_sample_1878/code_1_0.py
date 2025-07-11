import chess
import itertools

def find_checkmate_positions():
    """
    Calculates the number of legal chess positions where White has a King,
    Pawn, and Knight, and Black has only a King, and Black is in checkmate.
    """
    checkmate_count = 0

    # Generate all permutations of 4 unique squares out of 64 for the 4 pieces.
    # The order of pieces is: White King, White Knight, White Pawn, Black King.
    for squares in itertools.permutations(chess.SQUARES, 4):
        wk_square, wn_square, wp_square, bk_square = squares

        # --- Legality Check 1: Pawn Rank ---
        # A pawn cannot be on the 1st (rank 0) or 8th (rank 7) rank.
        if chess.square_rank(wp_square) in [0, 7]:
            continue

        # --- Legality Check 2: King Proximity ---
        # The two kings cannot be on adjacent squares.
        if chess.square_distance(wk_square, bk_square) <= 1:
            continue

        # --- Setup the board for further checks ---
        board = chess.Board(fen=None)  # Create an empty board
        board.set_piece_at(wk_square, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(wn_square, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_square, chess.Piece(chess.PAWN, chess.WHITE))
        board.set_piece_at(bk_square, chess.Piece(chess.KING, chess.BLACK))

        # --- Legality Check 3: Side not to move is not in check ---
        # For a position to be legal, the player whose turn it is NOT
        # cannot be in check. Since it's Black's turn in a checkmate,
        # we check if White is in check. If so, it's an illegal position.
        board.turn = chess.WHITE
        if board.is_check():
            continue

        # --- Final Check: Is it checkmate? ---
        # Set the turn to Black and check if the position is checkmate.
        board.turn = chess.BLACK
        if board.is_checkmate():
            checkmate_count += 1

    # The final count is the answer.
    print(f"The total number of constructible checkmate positions is: {checkmate_count}")

if __name__ == '__main__':
    find_checkmate_positions()