# The puzzle asks for the number of moves in which White can win, assuming optimal play.
#
# 1. White's first move: Be5
#    This move creates an unstoppable threat of moving the King to d4.
#    The d4 square is a winning square in King of the Hill.
#    The square d4 is not controlled by any of Black's pieces.
#
# 2. Black's move:
#    Black has no move that can stop White from playing Kd4 on the next turn.
#    Black cannot move a piece to attack d4 in one turn.
#    Let's assume Black plays a random move like ...Kb6.
#
# 3. White's second move: Kd4
#    White moves the King from e3 to d4. Since d4 is a central square, White wins the game.
#
# The sequence is: 1. Be5 ... 2. Kd4#
# White wins in 2 moves.

winning_moves = 2
print(winning_moves)