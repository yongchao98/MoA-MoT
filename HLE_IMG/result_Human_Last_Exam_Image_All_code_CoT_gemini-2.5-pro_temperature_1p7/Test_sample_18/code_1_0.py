def solve_chess_puzzle():
  """
  This script provides the solution to the chess puzzle.
  """
  hidden_piece = "White King"
  piece_location = "c3"
  mating_move = "1... Ra3#"

  print("The solution to the chess puzzle is as follows:")
  print(f"The hidden piece is the {hidden_piece}.")
  print(f"The hidden piece must be located on the square {piece_location}.")
  print(f"With the White King on {piece_location}, Black can mate in one move.")
  print(f"The mating move is: {mating_move}")

  print("\n--- Detailed Explanation ---")
  print("The White King is placed on c3, a square not attacked by any black piece.")
  print("Black plays 1... Ra3, moving the rook from a1 to a3.")
  print("This is checkmate because:")
  print("1. The King is in check from the Rook on a3.")
  print("2. The King has no escape squares:")
  print("   - Occupied squares: b2 (Queen), d2 (Pawn), b3 (Pawn), d3 (Bishop), c4 (Knight).")
  print("   - Attacked squares: b4 (by Rook on a3), d4 (by Rook on a3), c2 (by Queen on e2).")
  print("3. The check cannot be blocked or the Rook on a3 captured.")

solve_chess_puzzle()