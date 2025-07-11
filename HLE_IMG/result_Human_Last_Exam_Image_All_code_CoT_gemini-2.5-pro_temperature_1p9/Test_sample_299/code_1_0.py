def solve_chess_puzzle():
  """
  This function prints the solution to the mate-in-2 chess puzzle.
  """
  move1_white = "1. Qa7"
  move1_black = "   Kxf3"
  move2_white = "2. Qf2#"
  
  print("The checkmate in 2 moves is:")
  print(f"{move1_white}{move1_black}")
  print(f"{move2_white}")

solve_chess_puzzle()