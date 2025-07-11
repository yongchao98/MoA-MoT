def solve_chess_puzzle():
  """
  This function provides the solution to the chess puzzle.
  It identifies the minimum number of moves for black to force a checkmate
  and states the best move to initiate that sequence.
  """
  
  # The forced checkmate sequence is:
  # 1. ... R8a4+
  # 2. Ke4 Rf3#
  
  # The number of moves for black to deliver checkmate.
  number_of_moves = 2
  
  # The best starting move for black.
  best_move = "R8a4+"
  
  print(f"{number_of_moves}. Best move: {best_move}")

solve_chess_puzzle()