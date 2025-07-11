def solve_chess_puzzle():
  """
  This function prints the solution to the chess puzzle.
  The solution is based on retrograde analysis which proves that
  White's last move must have been d5-d6, allowing for a mate
  in one via an en passant capture.
  """
  move = "1... cxd6 e.p.#"
  print("The mating move is:")
  print(move)

solve_chess_puzzle()