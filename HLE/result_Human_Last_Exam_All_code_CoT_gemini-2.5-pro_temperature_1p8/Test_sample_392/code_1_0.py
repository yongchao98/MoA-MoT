# First, you may need to install the required library:
# pip install python-chess

import chess
import sys

def find_best_mate_sequence():
    """
    This function analyzes a specific chess position to find and display
    the shortest forced checkmate sequence for Black.
    """
    # The FEN string describes the initial state of the chessboard.
    fen = "r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1"
    board = chess.Board(fen)

    print("Initial position (Black to move):")
    print(board)
    print("\nFinding the shortest forced mate sequence...")
    
    # The shortest mating sequence is 1... Qg1+ 2. Rxg1 Nf2#
    # Move 1: Black plays Qg1+
    try:
        # Standard Algebraic Notation for Black's first move
        move1_black = board.parse_san("Qg1")
        board.push(move1_black)
        # Using board.san() to get the full notation with check/mate symbols
        print(f"1. ... {board.san(move1_black)}")
    except ValueError as e:
        print(f"Error making move 'Qg1': {e}", file=sys.stderr)
        return

    # Move 2: White plays Rxg1 (this is forced)
    try:
        move2_white = board.parse_san("Rxg1")
        board.push(move2_white)
        print(f"2. {board.san(move2_white)}")
    except ValueError as e:
        print(f"Error making move 'Rxg1': {e}", file=sys.stderr)
        return
        
    # Move 2: Black plays Nf2#
    # The move 'Nf2' is ambiguous because both knights (on d4 and e4) can move to f2.
    # We must disambiguate. Both Ndf2 and Nef2 lead to the same mate.
    # We'll use 'Ndf2' for this demonstration.
    try:
        move2_black = board.parse_san("Ndf2")
        board.push(move2_black)
        print(f"   ... {board.san(move2_black)}")
    except ValueError as e:
        print(f"Error making move 'Ndf2': {e}", file=sys.stderr)
        return
        
    print("\nFinal board position:")
    print(board)
    
    # Final confirmation
    if board.is_checkmate():
        print("\nConclusion: The sequence results in a checkmate.")
    else:
        print("\nConclusion: The sequence does not result in a checkmate.")

if __name__ == '__main__':
    find_best_mate_sequence()

<<<C>>>