def solve_heesch_number_puzzle():
  """
  This function determines and prints the Heesch number for the three
  identical polygons shown in the image.

  The analysis of the polygon's shape reveals that it is a space-filling tile.
  Its boundary is composed of three pairs of opposite sides. In each pair,
  the sides are complementary (one is a "bump", the other a "dent") and
  are parallel (related by a translation). Such a shape can tessellate the
  plane.

  By definition, the Heesch number of a shape that can tile the plane is infinity.
  """

  # The Heesch number for a plane-tiling polygon is infinity.
  # We use the Unicode character for the infinity symbol.
  heesch_number = "\u221E"

  # The problem asks for the Heesch number for each of the three identical polygons.
  answer_for_polygon_1 = heesch_number
  answer_for_polygon_2 = heesch_number
  answer_for_polygon_3 = heesch_number

  # We format the answer as a comma-separated string.
  final_answer = f"{answer_for_polygon_1}, {answer_for_polygon_2}, {answer_for_polygon_3}"
  
  print(final_answer)

solve_heesch_number_puzzle()