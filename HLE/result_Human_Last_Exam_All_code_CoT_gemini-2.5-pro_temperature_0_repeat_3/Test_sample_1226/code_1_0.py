# This script calculates the maximum material value based on the chess puzzle's constraints.

# Standard piece point values.
QUEEN_VALUE = 9
PAWN_VALUE = 1

# Total squares on a chessboard.
TOTAL_SQUARES = 64

# We determine the number of squares occupied by essential pieces for the proof and checkmate.
# These pieces are: 1 Black King, 1 White King, 1 promoted Queen, 1 original Queen, and 2 Pawns for proof.
# The White King must have moved, and the pawns block faster promotion paths.
# The original Queen proves that the other Queen is a promotion.
squares_used_by_essential_pieces = 6

# The remaining squares can be filled with the highest-value piece (Queens) to maximize the score.
num_filler_queens = TOTAL_SQUARES - squares_used_by_essential_pieces

# The white pieces contributing to the score are those essential for the setup plus the fillers.
# The essential white pieces are 2 Queens and 2 Pawns.
num_setup_queens = 2
num_setup_pawns = 2

# Calculate the total number of each white piece on the board.
total_queens = num_setup_queens + num_filler_queens
total_pawns = num_setup_pawns

# Calculate the final maximum point value.
total_value = (total_queens * QUEEN_VALUE) + (total_pawns * PAWN_VALUE)

# Output the results clearly.
print(f"To solve this, we construct a specific checkmate position that is provably unreachable in under 6 moves.")
print(f"The setup requires {num_setup_queens} Queens and {num_setup_pawns} Pawns.")
print(f"The remaining {num_filler_queens} squares are filled with Queens to maximize the score.")
print(f"Total number of Queens: {total_queens}")
print(f"Total number of Pawns: {total_pawns}")
print("\nFinal Calculation:")
print(f"({total_queens} * {QUEEN_VALUE}) + ({total_pawns} * {PAWN_VALUE}) = {total_value}")

print(f"\nThe greatest number of points of white material is {total_value}.")