def print_mating_sequence():
  """
  Prints the best sequence of moves for black to checkmate white, which is a smothered mate.
  The sequence is:
  1. Black moves Queen to g1, checking the White King.
  2. White is forced to capture the Queen with its Rook.
  3. Black moves Knight to f2, delivering a checkmate.
  """
  print("The best sequence of moves is a classic smothered mate:")
  print("1... Qg1+")
  print("2. Rxg1 Nf2#")

print_mating_sequence()