import chess
import itertools

def solve_chess_puzzle():
    """
    Calculates the number of possible legal checkmate positions where
    White has a King, Knight, and Pawn, and Black has only a King.

    Method:
    1. Generates every possible placement of the 4 pieces on the 64 squares.
    2. Filters out placements that are fundamentally illegal:
       - Pawns on the 1st or 8th rank.
       - Kings on adjacent squares.
    3. For each potentially legal placement, it sets up a board.
    4. It assumes it's Black's turn (as Black is the one checkmated).
    5. It checks if the position is a valid checkmate using the python-chess library.
    6. It counts and returns the total number of such positions.
    """
    checkmate_count = 0
    
    # Generate all unique placements of 4 pieces on a 64-square board.
    # The pieces are: Black King, White King, White Knight, White Pawn.
    piece_placements = itertools.permutations(range(64), 4)

    for placement in piece_placements:
        bk_sq, wk_sq, wn_sq, wp_sq = placement

        # Filter 1: A pawn cannot be on its back rank (rank 1 or 8).
        # In python-chess, ranks are 0-indexed (0 to 7).
        if chess.square_rank(wp_sq) in [0, 7]:
            continue
            
        # Filter 2: The kings can never be on adjacent squares.
        # We check the Chebyshev distance between the kings.
        if chess.square_distance(bk_sq, wk_sq) <= 1:
            continue

        # If the basic filters pass, set up the board for a full check.
        board = chess.Board(fen=None)  # Create an empty board
        board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))
        board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(wn_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_sq, chess.Piece(chess.PAWN, chess.WHITE))

        # In a checkmate position, it is the turn of the player who is mated.
        board.turn = chess.BLACK

        # Final Check: Is the position a checkmate?
        # The library function board.is_checkmate() efficiently verifies if
        # Black's king is in check and has no legal moves.
        if board.is_checkmate():
            checkmate_count += 1
            
    # The problem statement asked to output the final number.
    print(checkmate_count)

# Execute the main function to find and print the answer.
solve_chess_puzzle()
<<<16360>>>