import chess
import itertools

def find_checkmate_positions():
    """
    Calculates the number of legal checkmate positions with White (K, N, P)
    against Black (K).

    A position is a legal checkmate if:
    1. It is Black's turn to move.
    2. The Black king is in check.
    3. Black has no legal moves.
    4. The position is 'legal': the White king is not attacked by the Black king.
    5. White's pawn is not on the 1st or 8th rank.
    """
    count = 0
    squares = range(64)

    # Use itertools.permutations to efficiently generate all unique placements
    # of the 4 pieces on the 64 squares.
    # The tuple will contain square indices for (White King, Black King, White Knight, White Pawn)
    piece_placements = itertools.permutations(squares, 4)

    print("Analyzing all possible piece placements... (This may take a few minutes)")

    for placement in piece_placements:
        wk_sq, bk_sq, wn_sq, wp_sq = placement

        # Rule 5: A white pawn cannot be on the 1st or 8th rank.
        # Ranks are 0-indexed in the library (0=rank 1, 7=rank 8).
        wp_rank = chess.square_rank(wp_sq)
        if wp_rank == 0 or wp_rank == 7:
            continue

        # Set up the board with the current piece placement
        board = chess.Board.empty()
        board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))
        board.set_piece_at(wn_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
        board.set_piece_at(wp_sq, chess.Piece(chess.PAWN, chess.WHITE))

        # Rule 1: It must be Black's turn to move.
        board.turn = chess.BLACK

        # The python-chess library's functions check the remaining rules:
        # board.is_checkmate() handles rules 2 and 3.
        # board.is_valid() handles rule 4 by ensuring the side not to move (White) is not in check.
        if board.is_valid() and board.is_checkmate():
            count += 1

    # The final equation is the total sum of all found positions.
    print(f"Total number of legal checkmate positions = {count}")

if __name__ == '__main__':
    find_checkmate_positions()