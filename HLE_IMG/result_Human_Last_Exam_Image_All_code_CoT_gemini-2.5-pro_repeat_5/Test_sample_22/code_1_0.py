def find_unknots():
  """
  Analyzes the provided knot diagrams to identify which are isotopic to the unknot.

  Knot Analysis:
  - K_1: This is the 6_1 knot. It is not the unknot.
  - K_2: This is a complex diagram that can be simplified to the unknot using Reidemeister moves. It is the unknot.
  - K_3: This is a simple 2-crossing diagram of the unknot, undone by a Reidemeister II move. It is the unknot.
  - K_4: This is the figure-eight knot (4_1). It is not the unknot.
  - K_5: This is a known tricky diagram of the unknot. It is the unknot.
  - K_6: This is the cinquefoil knot (5_1). It is not the unknot.
  """

  # List of indices of knots that are isotopic to the unknot.
  unknot_indices = [2, 3, 5]

  print(unknot_indices)

find_unknots()