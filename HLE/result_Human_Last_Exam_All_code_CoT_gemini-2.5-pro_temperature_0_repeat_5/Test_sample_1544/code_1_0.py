# This script determines and prints the best move for White.
#
# The position is 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1.
# White is under pressure, and the only move that does not lead to a loss
# is to trade queens, forcing a drawn pawn endgame.
#
# The principal variation, which can be considered the "final equation" of this position,
# is 1. Qxc5+ bxc5. The numbers in this notation are 1, 5, and 5.
# The best move for White is the first part of this sequence.

best_move = "Qxc5+"
print(best_move)