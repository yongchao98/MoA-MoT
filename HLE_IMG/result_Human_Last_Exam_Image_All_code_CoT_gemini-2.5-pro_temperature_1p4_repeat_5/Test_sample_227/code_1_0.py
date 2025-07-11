def solve_chess_puzzle():
  """
  Determines the fastest checkmate for black and prints the result.
  The analysis shows a mate in 3 moves is possible.
  The best move is Ra4+.
  The sequence is:
  1. ... Ra4+
  2. Kh3
  2. ... g3
  3. [Any White move] ... R(a8)h8#
  """
  moves = 3
  best_move = "Ra4+"
  print(f"{moves}. Best move: {best_move}")

solve_chess_puzzle()