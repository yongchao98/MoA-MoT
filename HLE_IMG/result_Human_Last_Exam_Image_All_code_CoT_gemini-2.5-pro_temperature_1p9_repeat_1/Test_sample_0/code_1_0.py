def solve_chess_puzzle():
  """
  This function prints the solution to the chess puzzle.
  The puzzle is a mate in 2 for black, without moving the queens.

  The sequence is as follows:
  1. Nxc1: Black's knight at d1 captures White's bishop at c1. This removes the only piece
            that can defend the f2 square. White is forced to respond with a pawn move
            (either c4 or e6), as the king is trapped.
  2. Nf2#: Black's knight moves to f2, delivering checkmate. The White king cannot move
           as all its escape squares (g1, g2, h2) are controlled by the black queens,
           and the checking knight cannot be captured.
  """
  move1 = "Nxc1"
  move2 = "Nf2#"
  print(f"{move1}, {move2}")

solve_chess_puzzle()