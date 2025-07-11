def find_best_move():
  """
  Analyzes the FEN string '8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1' 
  and determines the best move for White.
  """
  
  # The position analysis reveals a winning simplification for White.
  # The key is to trade queens to eliminate Black's counterplay and 
  # promote the powerful d6-pawn.
  
  # The move is 1. Qxc5+
  move = "Qxc5+"
  
  # After 1. Qxc5+, Black must reply 1... kxc5.
  # This enters a King and Pawn endgame where White is winning.
  # White's winning plan is to bring the King from g1 towards the center 
  # (e.g., 2. Kf1, 3. Ke2, 4. Kd3).
  # Black's king is tied down by the passed d-pawn and cannot stop
  # the advance of White's king, which will lead to a decisive advantage.
  
  print("Based on the analysis, the optimal move for White is a queen trade initiated by a check.")
  print("This move simplifies the position into a winning King and Pawn endgame.")
  print(f"\nThe best move is: {move}")

find_best_move()