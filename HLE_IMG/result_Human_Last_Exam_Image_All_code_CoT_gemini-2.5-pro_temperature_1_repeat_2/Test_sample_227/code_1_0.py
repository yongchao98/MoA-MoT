def solve_chess_puzzle():
  """
  This function provides the solution to the given chess puzzle.
  It prints the minimum number of moves to force a checkmate and the best move for black.
  """
  # The analysis shows a forced checkmate in 2 moves.
  # 1... Ra4+
  # White's only legal reply is to block the check.
  # 2. Rd4
  # Black then captures the blocking rook for checkmate.
  # 2... Rxd4#
  number_of_moves = 2
  best_move = "Ra4+"
  print(f"{number_of_moves}. Best move: {best_move}")

solve_chess_puzzle()