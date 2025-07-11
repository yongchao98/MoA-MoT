def solve_chess_puzzle():
  """
  This function prints the solution to the chess mate-in-2 puzzle.
  """
  white_move_1 = "Ng6+"
  black_move_1 = "Rxg6"
  white_move_2_mate = "Rhg5#"

  print("The solution for the mate in 2 is:")
  # Print the first move for both sides
  print(f"1. {white_move_1} {black_move_1}")
  # Print the final checkmating move
  print(f"2. {white_move_2_mate}")

solve_chess_puzzle()