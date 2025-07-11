import chess
import itertools

def find_checkmate_positions():
    """
    Calculates the number of legal checkmate positions where White has a King,
    Pawn, and Knight, and Black has only a King.
    """
    checkmate_count = 0
    piece_types = ['WK', 'BK', 'WN', 'WP']
    total_permutations = 15249024 # 64*63*62*61
    processed_count = 0

    print("Starting analysis... This will take a few minutes.")

    # Iterate over all permutations of 4 unique squares out of 64
    for squares in itertools.permutations(range(64), len(piece_types)):
        wk_sq, bk_sq, wn_sq, wp_sq = squares
        
        processed_count += 1
        if processed_count % 1_000_000 == 0:
            print(f"Processed {processed_count:,} of {total_permutations:,} positions...")

        # --- Rule 1: Pawns cannot be on the 1st or 8th rank ---
        # A white pawn on rank 1 or 8 is illegal.
        # chess.square_rank(sq) gives rank index 0-7.
        wp_rank = chess.square_rank(wp_sq)
        if wp_rank == 0 or wp_rank == 7:
            continue

        # --- Rule 2: Kings cannot be adjacent ---
        # The distance between the squares must be greater than 1.
        if chess.square_distance(wk_sq, bk_sq) <= 1:
            continue

        # --- Setup the board for the current position ---
        board = chess.Board.empty()
        board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))
        board.set_piece_at(wn_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_sq, chess.Piece(chess.PAWN, chess.WHITE))

        # --- Rule 3: It must be Black's turn ---
        board.turn = chess.BLACK
        
        # --- Rule 4: The side not to move (White) cannot be in check ---
        # Since black only has a king, this is already covered by Rule 2.
        # We can proceed to check for mate.

        # --- Final Check: Is it checkmate? ---
        if board.is_checkmate():
            checkmate_count += 1

    print("\nAnalysis complete.")
    print("Total number of unique checkmate positions found:", checkmate_count)

if __name__ == '__main__':
    # Before running, make sure you have the python-chess library installed:
    # pip install python-chess
    find_checkmate_positions()
<<<4704>>>