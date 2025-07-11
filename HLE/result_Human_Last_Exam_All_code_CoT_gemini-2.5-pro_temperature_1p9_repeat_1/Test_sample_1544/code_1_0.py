# This script will print the best move for White in the given chess position.

# The FEN for the position is: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
# White's key advantages are the passed d-pawn and the exposed Black king.
# Black's main threat is the passed a-pawn.
# A direct attack on the king is the most effective strategy.
# The move Qf7+ initiates a forced sequence that leads to a win for White.
# For example: 1. Qf7+ Ka6 2. Qxd7 Qc1+ 3. Kh2 a4 4. Qc8+
# This line shows that White's attack combined with the d-pawn threat is overwhelming.

best_move = "Qf7+"
print(best_move)