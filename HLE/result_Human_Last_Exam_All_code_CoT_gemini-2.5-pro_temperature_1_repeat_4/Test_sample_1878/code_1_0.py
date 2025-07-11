import chess
import itertools

def find_checkmate_positions():
    """
    Calculates the number of legal checkmate positions with White (King, Pawn, Knight)
    vs. Black (King).

    A position is considered 'legal' if:
    1. The two kings are not on adjacent squares.
    2. The white pawn is not on the 1st or 8th rank.
    3. The side not to move (White) is not in check.

    A position is a checkmate if it's Black's turn, the Black King is in
    check, and Black has no legal moves.
    """
    checkmate_count = 0
    all_squares = chess.SQUARES

    # Generate all permutations of 4 unique squares for the 4 pieces.
    # The order is: White King, Black King, White Knight, White Pawn.
    # P(64, 4) = 15,249,024 total placements to check.
    piece_placements = itertools.permutations(all_squares, 4)

    for wk_sq, bk_sq, wn_sq, wp_sq in piece_placements:
        # --- Legality Pre-Check 1: Pawn Rank ---
        # A white pawn cannot be on the 1st or 8th rank.
        # Ranks are 0-indexed (0=rank 1, 7=rank 8).
        if chess.square_rank(wp_sq) in [0, 7]:
            continue

        # --- Legality Pre-Check 2: King Proximity ---
        # The two kings cannot be on adjacent squares.
        if chess.king_distance(wk_sq, bk_sq) <= 1:
            continue

        # --- Construct Board and Perform Final Checks ---
        board = chess.Board.empty()
        board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))
        board.set_piece_at(wn_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_sq, chess.Piece(chess.PAWN, chess.WHITE))

        # It must be Black's turn in a checkmate position against Black.
        board.turn = chess.BLACK

        # --- Legality Check 3 & Checkmate Verification ---
        # board.is_valid() checks for several conditions, including that the
        # side NOT to move (White) is not in check.
        # board.is_checkmate() checks if the side TO move (Black) is checkmated.
        # We need both conditions to be true for a valid checkmate.
        if board.is_valid() and board.is_checkmate():
            checkmate_count += 1

    print(f"The total number of constructible checkmate positions is: {checkmate_count}")

if __name__ == '__main__':
    # This calculation can take a few minutes to run.
    find_checkmate_positions()