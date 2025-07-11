import sys

def solve_tropical_moduli_question():
  """
  This function provides the answer to the theoretical question about tropical moduli spaces.
  The reasoning is provided in the text explanation.
  """
  # Part (a): Minimum number of vertices for non-empty M_trop(g,A).
  # As explained, if the space is non-empty, a stable graph with one vertex exists.
  part_a = "1"

  # Part (b): Is M_trop(0,A) always a simplicial fan?
  # Yes, this is a standard result in tropical geometry.
  part_b = "yes"

  # Part (c): For g > 0, is M_trop(g,A) a tropical variety? If not, is it a polyhedral complex?
  # No, because for g >= 2 it is not pure-dimensional.
  # Yes, it is a polyhedral complex by construction.
  part_c = "no, yes"

  # Formatting the final answer string
  final_answer = f"(a) {part_a}; (b) {part_b}; (c) {part_c}"

  print(final_answer)

solve_tropical_moduli_question()