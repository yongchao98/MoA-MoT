# This script prints the solution to the chess puzzle.
# The puzzle is a mate in 2 for Black.
#
# The sequence is determined by the following logic:
# 1. Black's first move is Knight to g3, which is check (Ng3+).
# 2. This move forces White's King to h1, as all other escape squares are
#    blocked or attacked, and capturing the Knight is an illegal move.
# 3. Black's second move is Rook to h8, which is checkmate (Rh8#).
#    The king has no escape squares, the check cannot be blocked, and the
#    Rook cannot be captured.
#
# The final sequence of black moves is printed below.

move1_number = 3
move2_number = 8
move1 = f"Ng{move1_number}+"
move2 = f"Rh{move2_number}#"

print("The sequence is:")
print(f"{move1}, {move2}")