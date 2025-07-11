def find_unknots():
  """
  Analyzes the six knots K_1 to K_6 to identify which are isotopic to the unknot.

  The analysis is as follows:
  - K_1 is the 6_2 knot, not an unknot.
  - K_2 is a complicated diagram of the unknot, simplifiable by Reidemeister II moves.
  - K_3 can be simplified by shrinking one of its loops to a point.
  - K_4 is the figure-eight knot, not an unknot.
  - K_5 is isotopic to the figure-eight knot via a Reidemeister I move.
  - K_6 is the square knot (a composite knot), not an unknot.

  Therefore, only knots K_2 and K_3 are unknots.
  """
  unknot_indices = [2, 3]
  print(unknot_indices)

find_unknots()