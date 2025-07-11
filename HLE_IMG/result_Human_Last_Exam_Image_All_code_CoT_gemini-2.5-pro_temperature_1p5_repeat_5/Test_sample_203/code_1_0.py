def solve_heesch_numbers():
  """
  Calculates and prints the Heesch numbers for the three polygons in the image.
  A Heesch number of infinity means the polygon tiles the plane.
  """
  
  # Heesch number for the first (top) polygon.
  # It is centrally symmetric, so it tiles the plane.
  heesch_num_1 = "∞"
  
  # Heesch number for the second (middle) polygon.
  # Its boundary tabs and notches are not complementary, so it cannot form a full layer.
  heesch_num_2 = "0"

  # Heesch number for the third (bottom) polygon.
  # Its opposite sides are complementary, so it tiles the plane.
  heesch_num_3 = "∞"
  
  # The problem asks for the answers separated by commas.
  print(f"{heesch_num_1},{heesch_num_2},{heesch_num_3}")

solve_heesch_numbers()