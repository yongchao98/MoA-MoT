def solve_chess_puzzle():
  """
  This function provides the solution to the chess puzzle.
  It identifies the minimum number of moves for black to force a checkmate
  and the best move to achieve it.
  """
  # The number of moves for the forced checkmate.
  number_of_moves = 3
  
  # The best move for black in Standard Algebraic Notation.
  best_move = "R3a4+"
  
  print(f"{number_of_moves}. Best move: {best_move}")

solve_chess_puzzle()