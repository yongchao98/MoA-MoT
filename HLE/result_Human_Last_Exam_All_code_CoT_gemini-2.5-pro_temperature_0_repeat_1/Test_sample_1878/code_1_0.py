import chess
import itertools

def count_kpn_vs_k_checkmates():
    """
    Calculates the number of legal checkmate positions where White has a King,
    Pawn, and Knight, and Black has only a King.

    A position is considered 'legal' if:
    1. No two pieces occupy the same square.
    2. The pawn is not on the 1st or 8th rank.
    3. The two kings are not on adjacent squares.

    This script requires the 'python-chess' library.
    You can install it using: pip install python-chess
    """
    checkmate_count = 0
    all_squares = range(64)

    # Generate all permutations of 4 unique squares for the 4 pieces.
    # The order of pieces is: White King, White Knight, White Pawn, Black King.
    total_permutations = 64 * 63 * 62 * 61 # Just for info, not used in loop
    
    # Using itertools.permutations is an efficient way to handle this.
    for squares_permutation in itertools.permutations(all_squares, 4):
        wk_square, wn_square, wp_square, bk_square = squares_permutation

        # --- Filter 1: Pawn Rank Rule ---
        # A pawn cannot be on its army's first rank or the promotion rank.
        pawn_rank = chess.square_rank(wp_square)
        if pawn_rank == 0 or pawn_rank == 7:
            continue

        # --- Filter 2: King Proximity Rule ---
        # The two kings cannot be on adjacent squares.
        if chess.square_distance(wk_square, bk_square) <= 1:
            continue

        # --- Construct the Board and Check for Mate ---
        # Create an empty board
        board = chess.Board(fen=None)

        # Place the pieces on the board
        board.set_piece_at(wk_square, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(wn_square, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_square, chess.Piece(chess.PAWN, chess.WHITE))
        board.set_piece_at(bk_square, chess.Piece(chess.KING, chess.BLACK))

        # To check for checkmate against Black, it must be Black's turn.
        board.turn = chess.BLACK

        # The is_checkmate() function verifies that Black is in check
        # and has no legal moves. This also implicitly confirms that the
        # position is legal from the current turn's perspective.
        if board.is_checkmate():
            checkmate_count += 1
            
    # The final equation is simply the total count.
    print(f"Total checkmate positions: {checkmate_count}")

if __name__ == '__main__':
    # Note: This script is computationally intensive and may take a few minutes to run.
    count_kpn_vs_k_checkmates()