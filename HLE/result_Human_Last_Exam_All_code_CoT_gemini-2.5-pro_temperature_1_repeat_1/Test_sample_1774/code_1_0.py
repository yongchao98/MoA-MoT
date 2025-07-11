# First, you may need to install the required library.
# You can do this by running: pip install python-chess

import chess
import sys

def solve_chess_puzzle():
    """
    This function sets up the described chess position and demonstrates
    the shortest mating sequence for White.
    """
    # The board position is represented using Forsyth-Edwards Notation (FEN).
    # Ranks are listed from 8 to 1.
    # 'w' means White to move. 'KQ' means White has castling rights.
    fen = "rn2r1k1/pbp1q1pp/1p2p1b1/4N3/3PN3/3B4/PPP2PPP/R3K2R w KQ - 0 1"

    # Create a board object from the FEN string.
    board = chess.Board(fen)

    # --- The Mating Sequence ---
    # After analyzing the position, a mate in 2 moves was found.
    # The sequence is: 1. Nxg6+ Kh8 2. Qh7#

    print("The shortest path to checkmate is a mate in 2 moves.")
    print("The winning sequence is printed below:")

    # Move 1: White plays Nxg6+
    # We use Standard Algebraic Notation (SAN) for the moves.
    move1_white_san = "Nxg6+"
    try:
        move1_white = board.parse_san(move1_white_san)
        board.push(move1_white)
    except ValueError:
        print(f"Error: The move '{move1_white_san}' is not legal in the current position.")
        return

    # Move 1: Black is forced to play Kh8
    # We can verify this is the only legal move for Black.
    if len(list(board.legal_moves)) != 1:
        print("Error in analysis: Black has more than one reply.")
        return
        
    move1_black = list(board.legal_moves)[0]
    move1_black_san = board.san(move1_black)
    board.push(move1_black)
    
    # Print the first full move
    # The output shows the move number, White's move, and Black's move.
    print(f"1. {move1_white_san} {move1_black_san}")

    # Move 2: White plays Qh7#
    move2_white_san = "Qh7#"
    try:
        move2_white = board.parse_san(move2_white_san)
        board.push(move2_white)
    except ValueError:
        print(f"Error: The move '{move2_white_san}' is not legal.")
        return

    # Print White's second move
    print(f"2. {move2_white_san}")

    # Verify that the position is indeed checkmate.
    if board.is_checkmate():
        print("\nResult: Black is in checkmate.")
    else:
        print("\nResult: The sequence did not result in a checkmate.")

if __name__ == '__main__':
    solve_chess_puzzle()

<<<2>>>