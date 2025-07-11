def solve_chess_puzzle():
  """
  This function prints the solution to the chess puzzle.
  The puzzle is a mate in 2 for black, without moving the queens.
  The sequence is determined by analyzing the possible moves.

  Move 1: Rf2+
  This forces the white king to h1, as the rook is protected by the knight on d1,
  and the f1 square is controlled by the black queen on b1.

  Move 2: Rf1#
  This is checkmate. The king on h1 cannot move to g1 (queen on b1),
  g2 (bishop on h3), or h2 (pawn). The rook on f1 is protected by the queen on b1.
  """
  # The sequence notation includes the numbers for the squares (ranks and files).
  # The first move is Rook to f2, check.
  move1 = "Rf2+"
  # The second move is Rook to f1, checkmate.
  move2 = "Rf1#"

  print(f"{move1}, {move2}")

solve_chess_puzzle()