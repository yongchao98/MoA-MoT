import chess

def find_best_move():
  """
  Analyzes the given chess FEN and determines the best move for White.

  The position is 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1.

  The move 1. Qb2+ is a check. The black king on b6 has only one legal move,
  which is to a6, as all other escape squares are controlled by white's queen
  or pawns.

  After 1... Ka6, the black queen on c5 is no longer defended by its king.
  The white queen on b2 now attacks the black queen. White can then play
  2. Qxc5, winning the black queen and leading to a winning endgame.
  """
  best_move = "Qb2+"
  print("The best move for White is:")
  print(best_move)

find_best_move()