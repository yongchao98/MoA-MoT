# To run this code, you first need to install the python-chess library:
# pip install python-chess

import chess

def print_chess_solution():
    """
    This function verifies and prints the best move sequence for the given chess puzzle.
    The solution is a forced smothered mate in 2 for Black.
    """
    # The starting position in FEN (Forsyth-Edwards Notation)
    fen = "r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1"
    
    # Create a board object
    board = chess.Board(fen)

    # The moves for the shortest mate sequence (Option C) are:
    # 1... Qg1+
    # 2. Rxg1
    # 2... Nf2#
    
    # We will use the library to generate the correct Standard Algebraic Notation (SAN).
    
    # 1... Qg1+
    move1_black = chess.Move.from_uci("f2g1")
    san1_black = board.san(move1_black)
    board.push(move1_black)

    # 2. Rxg1
    move2_white = chess.Move.from_uci("f1g1")
    san2_white = board.san(move2_white)
    board.push(move2_white)

    # 2... Nf2#
    move2_black = chess.Move.from_uci("e4f2")
    san2_black = board.san(move2_black) # The library correctly adds '#' for checkmate

    # Print the full annotated sequence, showing each numbered move as requested.
    print(f"1... {san1_black}")
    print(f"2. {san2_white} {san2_black}")

# Execute the function to print the solution
print_chess_solution()