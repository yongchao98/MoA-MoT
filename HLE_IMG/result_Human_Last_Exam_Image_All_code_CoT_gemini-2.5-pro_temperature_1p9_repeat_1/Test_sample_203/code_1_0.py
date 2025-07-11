def solve_heesch_number():
  """
  This function determines and prints the Heesch number for the given polygons.

  The polygon in the image is a centrally-symmetric, modified hexagon.
  Any polygon with central symmetry can tile the plane.
  A shape that can tile the plane has a Heesch number of infinity.
  The image shows three identical polygons.
  Therefore, the Heesch number for each is infinity.
  """
  
  # Unicode for infinity symbol
  infinity_symbol = u'\u221e'
  
  # The Heesch numbers for the three identical polygons
  heesch_numbers = [infinity_symbol, infinity_symbol, infinity_symbol]
  
  # Print the result separated by commas
  print(','.join(heesch_numbers))

solve_heesch_number()