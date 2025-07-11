import chess
import itertools

def count_constructible_checkmates():
    """
    Calculates the number of legal chess positions where White has a King,
    Pawn, and Knight, Black has only a King, and Black is checkmated.
    """
    checkmate_count = 0
    total_permutations = 64 * 63 * 62 * 61
    
    # We iterate through all permutations of 4 unique squares out of 64.
    # Each permutation represents a potential placement for the four pieces.
    print(f"Analyzing all {total_permutations:,} possible piece placements. This will take several minutes...")
    
    # Using itertools.permutations is an efficient way to generate all unique placements
    for p in itertools.permutations(range(64), 4):
        # Assign squares from the permutation to the pieces
        wk_square, wn_square, wp_square, bk_square = p

        # --- Legality Check 1: White Pawn Rank ---
        # A white pawn cannot be on the 1st rank (rank index 0) or the 8th rank (rank index 7).
        # chess.square_rank() returns 0 for the 1st rank and 7 for the 8th.
        if chess.square_rank(wp_square) in (0, 7):
            continue

        # --- Legality Check 2: King Proximity ---
        # The two kings cannot be on adjacent squares. The distance must be greater than 1.
        # chess.square_distance measures the Chebyshev distance (max of file/rank differences).
        if chess.square_distance(wk_square, bk_square) <= 1:
            continue

        # --- Setup the Board ---
        # If the basic placement is plausible, create a board object to check full chess rules.
        board = chess.Board(fen=None)  # Start with an empty board
        board.set_piece_at(wk_square, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(wn_square, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_square, chess.Piece(chess.PAWN, chess.WHITE))
        board.set_piece_at(bk_square, chess.Piece(chess.KING, chess.BLACK))

        # --- Legality Check 3: Side to Move & Position Validity ---
        # For a checkmate to be registered, it must be the checkmated player's turn.
        board.turn = chess.BLACK
        
        # The 'is_valid()' method checks core chess rules. Crucially, it ensures
        # the side NOT to move (White) is not in check. A position where both
        # kings are in check is illegal and could not arise in a real game.
        if not board.is_valid():
            continue

        # --- Final Condition: Is it Checkmate? ---
        # If the position is fully legal, we finally check if Black is in checkmate.
        if board.is_checkmate():
            checkmate_count += 1
            
    print("Calculation complete.")
    # The prompt requests that the final output line prints the resulting number.
    print(f"The number of constructible checkmates is: {checkmate_count}")

# Execute the function
count_constructible_checkmates()