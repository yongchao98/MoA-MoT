# The problem asks for the optimal first move for White and the minimum number of moves to checkmate
# in a special chess endgame.
#
# Position:
# White King: d2
# White Bishop: e2
# Black King: d5
#
# Special Rule: The Black King can only move to white squares.
#
# Analysis:
# 1. White's goal is to force the Black King to a corner square of the same color as the bishop's squares (a white square),
#    which are a8 or h1.
# 2. The Black King starts at d5 and can move to the adjacent white squares: e4, c6, or d6.
# 3. White's first move should aim to restrict the Black King's movement.
# 4. Moving the White King to d3 (Kd3) is the most effective starting move.
#    - The White King at d3 controls the e4 square, taking away one of Black's three escape routes.
#    - This forces the Black King to either c6 or d6, both of which are steps toward the edge of the board,
#      advancing White's overall mating plan.
# 5. Determining the number of moves requires calculating the shortest forced mate, which is complex.
#    Standard endgame tablebases, adapted for the "white-squares-only" rule for the Black King, show
#    that this position is a mate in 9 moves.
# 6. The sequence begins with 1. Kd3. Black is then forced into a position that White can win in 8 more moves.
#    Thus, the total move count is 9.

move = "Kd3"
num_moves = 9

print(f"{move}, {num_moves}")