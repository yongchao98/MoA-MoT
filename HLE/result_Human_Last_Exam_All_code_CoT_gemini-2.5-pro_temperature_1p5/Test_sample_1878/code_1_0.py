import chess
import itertools

def count_checkmate_positions():
    """
    Calculates the number of legal checkmate positions with White (King, Knight, Pawn)
    vs Black (King).
    
    A position is considered legal if:
    1. Pieces are on distinct squares.
    2. The White pawn is not on rank 1 or 8.
    3. The two kings are not adjacent.
    4. White's king is not in check (implicit in this setup).

    A position is a checkmate if, with Black to move, Black's king is in
    check and Black has no legal moves.
    """
    checkmate_count = 0
    all_squares = range(64)

    # 1. Generate unique combinations of 4 squares for the 4 pieces.
    for squares_combo in itertools.combinations(all_squares, 4):
        # 2. For each combination, generate all piece assignments.
        for p in itertools.permutations(squares_combo):
            wk_sq, bk_sq, wn_sq, wp_sq = p

            # --- Legality Pre-checks ---

            # A) The White Pawn cannot be on the first or last rank (ranks 1 or 8).
            #    Ranks are 0-indexed in the library, so we check for rank 0 or 7.
            if chess.square_rank(wp_sq) in [0, 7]:
                continue

            # B) The Kings cannot be on adjacent squares.
            if chess.square_distance(wk_sq, bk_sq) <= 1:
                continue

            # --- Board Setup and Checkmate Verification ---

            # Create an empty board.
            board = chess.Board(fen=None)
            
            # Place the pieces on the board.
            board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
            board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))
            board.set_piece_at(wn_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
            board.set_piece_at(wp_sq, chess.Piece(chess.PAWN, chess.WHITE))
            
            # Set the turn to Black to check if it's checkmate.
            board.turn = chess.BLACK
            
            # C) Check if the position is checkmate.
            if board.is_checkmate():
                checkmate_count += 1
                
    # Final Output
    print(f"The total number of legal checkmate positions is: {checkmate_count}")

if __name__ == '__main__':
    count_checkmate_positions()