def solve_chess_puzzle():
  """
  This function provides the solution to the chess puzzle.
  The reasoning is as follows:
  1. Black's best move is 1... Rf3+, putting the White king in check.
  2. White's only legal response is 2. Kg5. All other squares are attacked.
  3. Black then plays 2... Rg3#, which is a checkmate.
  Therefore, it is a forced checkmate in 2 moves for Black.
  """
  number_of_moves = 2
  best_move = "Rf3+"
  print(f"{number_of_moves}. Best move: {best_move}")

solve_chess_puzzle()