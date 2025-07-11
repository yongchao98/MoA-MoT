def solve_chess_puzzle():
  """
  This function provides the solution to the mate-in-2 chess puzzle for Black.
  It identifies the sequence of moves that lead to a checkmate regardless of White's response,
  while adhering to the constraint of not moving the black queens.
  """
  
  # Black's first move: Knight to d2, delivering check.
  move1 = "Nd2+"
  
  # After White's only legal response (Bf2), Black's second move:
  # Rook captures f2, delivering checkmate.
  move2 = "Rxf2#"
  
  print("The sequence of moves for Black is:")
  print(f"1. {move1}")
  print(f"2. {move2}")

solve_chess_puzzle()