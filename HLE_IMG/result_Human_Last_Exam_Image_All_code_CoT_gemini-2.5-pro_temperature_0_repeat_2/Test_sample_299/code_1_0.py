def solve_chess_mate_in_two():
  """
  This function prints the solution for the given chess puzzle.
  The puzzle is a checkmate in 2 moves for White.
  """
  
  # White's first move, a check
  white_move_1 = "Ng6+"
  
  # Black's forced response
  black_move_1 = "Kd3"
  
  # White's second move, delivering checkmate
  white_move_2 = "Qd1#"
  
  print("The solution for the mate in 2 is:")
  # The problem asks to output each number in the final equation.
  # I will print the move numbers and the moves in algebraic notation.
  print(f"1. {white_move_1} {black_move_1}")
  print(f"2. {white_move_2}")

solve_chess_mate_in_two()