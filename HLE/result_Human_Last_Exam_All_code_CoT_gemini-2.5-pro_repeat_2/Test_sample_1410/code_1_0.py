def solve_chess_puzzle():
  """
  Prints the solution to the chess puzzle.
  The solution is a mate in 3 for White.
  """
  moves = 3
  winning_sequence = "1. Qe6+ Kh8 2. Qxe8+ Rxe8 3. Rxe8#"

  print(f"White can win in {moves} moves.")
  print("The winning sequence against optimal defense is:")
  print(winning_sequence)

solve_chess_puzzle()