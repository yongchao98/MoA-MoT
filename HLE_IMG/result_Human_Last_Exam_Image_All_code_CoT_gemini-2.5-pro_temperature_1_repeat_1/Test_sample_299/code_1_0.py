def solve_chess_mate_in_two():
  """
  This function prints the solution for the given chess puzzle.
  The solution is a 2-move checkmate for White.
  """
  
  # The move sequence for the checkmate
  white_move_1 = "1."
  white_piece_move_1 = "Qa1+"
  
  black_move_1 = "Kxf3"
  
  white_move_2 = "2."
  white_piece_move_2 = "Ng6#"
  
  print("The checkmate in 2 moves is:")
  print(f"{white_move_1} {white_piece_move_1} {black_move_1}")
  print(f"{white_move_2} {white_piece_move_2}")

solve_chess_mate_in_two()