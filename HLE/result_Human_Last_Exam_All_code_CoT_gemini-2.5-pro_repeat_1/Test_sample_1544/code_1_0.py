# This script will print the best move for White in the given chess position.
# The FEN for the position is: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
#
# Analysis Summary:
# 1. White's main advantage is the passed pawn on d6.
# 2. Black's queen on c5 is active and prevents the pawn from promoting.
# 3. The winning move for White is Qe5.
# 4. Qe5 pins the black queen to the black king, forcing a king move.
# 5. After the king moves (e.g., to b7), White plays Qxc5.
# 6. After the forced recapture Kxc5, White plays d7, and the pawn will promote to a new queen, leading to a win.

best_move = "Qe5"
print(f"The best move for White is {best_move}.")
