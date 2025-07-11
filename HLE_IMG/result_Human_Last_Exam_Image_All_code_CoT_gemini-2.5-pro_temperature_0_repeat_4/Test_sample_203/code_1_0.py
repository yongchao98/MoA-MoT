import math

def get_heesch_numbers():
  """
  This function determines the Heesch number for the polygons in the image.
  
  The analysis of the polygon's shape reveals that it is centrally symmetric.
  A key theorem in geometry states that any centrally symmetric polygon can tile the plane.
  The Heesch number for any shape that tiles the plane is infinity.
  The image shows three identical copies of this polygon.
  Therefore, the Heesch number for each of the three polygons is infinity.
  """
  
  # Using the infinity symbol '∞' as specified in the problem description.
  heesch_number_1 = "∞"
  heesch_number_2 = "∞"
  heesch_number_3 = "∞"
  
  # The problem asks for the answers in order, separated by commas.
  print(f"{heesch_number_1}, {heesch_number_2}, {heesch_number_3}")

get_heesch_numbers()