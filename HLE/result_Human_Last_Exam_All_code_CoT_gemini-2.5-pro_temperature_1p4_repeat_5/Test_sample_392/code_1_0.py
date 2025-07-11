def print_chess_solution():
  """
  Prints the correct, forced checkmating sequence for the given chess puzzle.
  """
  # The best sequence is a forced Smothered Mate as described in Option E.
  # We format it using Standard Algebraic Notation.
  solution_sequence = [
      "1. ... Nf2+",
      "2. Kg1 Nh3+",
      "3. Kh1 Qg1+",
      "4. Rxg1 Nf2#"
  ]
  
  print("The best sequence of moves for black and white that ends up with checkmating white is:")
  for move in solution_sequence:
    print(move)

print_chess_solution()
<<<E>>>