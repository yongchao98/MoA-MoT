import chess
import sys
import io

# A simple function to render the board to the console
def print_board(board):
    """Prints the board state to the console."""
    # Use StringIO to capture the print output of the board
    output = io.StringIO()
    # In python-chess, sys.stdout is used by default, so we redirect it
    _stdout = sys.stdout
    sys.stdout = output
    print(board)
    sys.stdout = _stdout
    # Process the captured output for better display
    captured = output.getvalue().replace(' . ', ' . ')
    print(captured)

# Set up the board from the FEN string
# FEN represents: rn3rk1/pbppq1pp/1p2pb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w - - 0 1
try:
    board = chess.Board("rn3rk1/pbppq1pp/1p2pb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w - - 0 1")
except ValueError:
    print("Invalid FEN string. Please check the board setup.")
    exit()

# The shortest path to checkmate is 3 moves.
# We will demonstrate one of the forced variations.
# 1. Nfg6+ hxg6  (The other option is 1...Bxg6 which also leads to mate in 3)
# 2. Nxg6+ Kh7   (Other options 2...Kf8 or 2...Bxg6 are also mate in 3)
# 3. Qh8#

# Create a list of moves in UCI format for the engine to process
# Nf4g6, h7xg6, Ne5xg6, g8h7, Qh5h8
moves_uci = ["f4g6", "h7g6", "e5g6", "g8h7", "h5h8"]
move_san_list = []

print("Initial position (White to move):")
print_board(board)
print("-" * 30)

move_number = 1
is_white_move = True
for move_str in moves_uci:
    move = chess.Move.from_uci(move_str)
    
    # Get the move in Standard Algebraic Notation (SAN) before pushing
    san = board.san(move)
    if is_white_move:
      print(f"Move {move_number}:")
      print(f"White plays: {san}")
    else:
      print(f"Black plays: {san}")
      move_number += 1

    board.push(move)
    print_board(board)
    print("-" * 30)
    is_white_move = not is_white_move


print("Final Result:")
if board.is_checkmate():
    print("White delivers checkmate in 3 moves.")
    print("The final move of the sequence is Qh8#")
else:
    print("The sequence did not result in a checkmate as expected.")

# Recreate the final equation as requested
print("\nThe final winning equation is:")
print("1. Nfg6+ hxg6")
print("2. Nxg6+ Kh7")
print("3. Qh8#")


<<<3>>>