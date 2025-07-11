def solve_heesch_number():
  """
  Determines the Heesch number for the polygons in the image.

  The image shows three identical polygons. The background pattern is a tessellation
  (a tiling of the plane) made from this same polygon.

  According to the definition of the Heesch number, if a shape can tile the plane,
  its Heesch number is considered to be infinity.

  Since the polygon shown tiles the plane, the Heesch number for each of the three
  identical polygons is infinity.
  """
  
  # Heesch number for the first polygon
  heesch_1 = "∞"
  
  # Heesch number for the second polygon
  heesch_2 = "∞"
  
  # Heesch number for the third polygon
  heesch_3 = "∞"
  
  # Print the answers in order, separated by commas.
  print(f"{heesch_1}, {heesch_2}, {heesch_3}")

solve_heesch_number()