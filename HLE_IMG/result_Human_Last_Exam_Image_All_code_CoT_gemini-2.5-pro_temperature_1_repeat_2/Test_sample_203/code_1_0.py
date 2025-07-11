def find_heesch_numbers():
  """
  This function determines and prints the Heesch number for the three
  polygons shown in the figure.

  The analysis of the polygon's shape reveals that it is a modified hexagon
  designed to tile the plane. It has pairs of opposite sides with
  complementary bumps and dents, allowing them to interlock perfectly.

  A polygon that can tile the plane can be surrounded by an infinite number of
  layers of itself. By definition, its Heesch number is infinity.

  Since all three polygons in the image are identical, they all have the
  same Heesch number.
  """
  
  # Heesch number for a single polygon that tiles the plane.
  # We use a string to represent infinity as requested by the problem's notation.
  heesch_number = "infinity"
  
  # The problem asks for the Heesch number of each of the three polygons,
  # separated by commas.
  polygon1_heesch = heesch_number
  polygon2_heesch = heesch_number
  polygon3_heesch = heesch_number
  
  print(f"{polygon1_heesch}, {polygon2_heesch}, {polygon3_heesch}")

find_heesch_numbers()