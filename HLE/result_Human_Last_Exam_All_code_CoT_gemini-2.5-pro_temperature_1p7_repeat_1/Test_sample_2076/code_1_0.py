# The goal is to find the minimum number of pieces for a Diagonal Corridor Mate
# with the White King on a1 and the Black King on h8, delivered by White.
# If multiple solutions exist with the same number of pieces, the one with the
# minimum total piece value is chosen.

# 1. The Black King on h8 has three escape squares: g8, h7, and g7.
# 2. A White piece must deliver check along the a1-h8 diagonal. This can be a
#    Bishop or a Queen.
# 3. We must control all three escape squares.

# Solution Analysis:
# - A 3-piece solution is the minimum possible. A checkmate requires at least two
#   pieces (one to check, one to control an escape square), and since the kings
#   cannot guard each other's escape squares, a third piece is needed to control
#   the remaining escape routes.
# - Two 3-piece solutions were identified:
#   1) White Queen, Black Pawn, Black Pawn. The Queen checks h8 and controls g7.
#      The two Black Pawns block g8 and h7.
#      Value = Queen(9) + Pawn(1) + Pawn(1) = 11.
#   2) White Bishop, White Knight, Black Pawn. The Bishop checks h8. The White
#      Knight (e.g., on f5) controls g7 and h7. The Black Pawn blocks g8.
#      Value = Bishop(3) + Knight(3) + Pawn(1) = 7.

# The second solution has the minimum total piece value (7).

# The pieces are:
# White: Bishop, Knight
# Black: Pawn

# Formatting the output: White pieces first, sorted by value, then Black pieces.
# Since Bishop and Knight have the same value (3), we list them alphabetically.
white_pieces = ["White Bishop", "White Knight"]
black_pieces = ["Black Pawn"]

# Combine the lists for the final output string
final_pieces_list = white_pieces + black_pieces
final_answer = ", ".join(final_pieces_list)

print(final_answer)