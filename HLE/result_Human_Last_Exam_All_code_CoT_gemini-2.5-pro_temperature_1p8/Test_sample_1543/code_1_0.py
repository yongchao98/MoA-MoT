# This problem is a "mate in 2" chess puzzle.
# White's first move is not a check, but it creates a situation
# where, no matter what Black plays, White can deliver checkmate on the second move.

# White's first move: 1. Qa6
# After this move, for any of Black's possible responses, White has a mating move.
# Therefore, White checkmates on the second move.

moves_to_win = 2

print("The minimal number of moves for White to win is:")
print(moves_to_win)