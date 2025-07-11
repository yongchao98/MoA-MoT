# The problem asks for the minimal number of moves by White to force a win,
# assuming optimal play from both sides.
# Our step-by-step analysis reveals a forced mate in 3.

# 1. Ag4 (Threatening Qd8#)
#    Black's best defense is Cd7. Any other move leads to mate in 2.
# 2. Qf5 (Now threatening Qf8#)
#    Black's only defense is Cf8.
# 3. Qxf8# (Checkmate)

# The number of moves for White in this sequence is 3.
minimal_moves_to_win = 3

print(minimal_moves_to_win)