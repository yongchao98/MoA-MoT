# The chess position is defined by the FEN: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
# White's best move is to trade queens, simplifying into a winning King and Pawn endgame.
# This move forces Black's king to recapture, leaving it out of position to stop White's other strategic threats.
# The move is Queen from f2 captures the piece on c5, delivering a check.
# In long algebraic notation, this is written as Qf2xc5+.
# This notation specifies the starting square (f2) of the piece moving, making it unambiguous.

move = "Qf2xc5+"
start_rank = 2
end_rank = 5

# Printing the components of the move as requested.
# 'Q' for Queen, 'f' for the starting file, '2' for the starting rank, 'x' for capture,
# 'c' for the destination file, '5' for the destination rank, '+' for check.
print(f"The optimal move for White is {move}.")
print(f"This involves moving the Queen from rank {start_rank} to rank {end_rank}.")
