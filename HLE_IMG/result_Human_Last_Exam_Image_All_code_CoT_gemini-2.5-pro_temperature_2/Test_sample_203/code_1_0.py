import math

def get_heesch_numbers():
  """
  Calculates and returns the Heesch numbers for the polygons in the image.

  Based on geometric analysis, the polygon is centrally symmetric.
  Any centrally symmetric polygon can tile the plane.
  If a polygon tiles the plane, its Heesch number is defined as infinity.
  The image shows three identical polygons.
  """
  
  # For a polygon that tiles the plane, the Heesch number is infinity.
  heesch_number_1 = math.inf
  heesch_number_2 = math.inf
  heesch_number_3 = math.inf
  
  # The problem asks for the answer as a comma-separated list.
  # We represent infinity with the string "infinity".
  answers = ["infinity", "infinity", "infinity"]
  
  print(", ".join(answers))

get_heesch_numbers()