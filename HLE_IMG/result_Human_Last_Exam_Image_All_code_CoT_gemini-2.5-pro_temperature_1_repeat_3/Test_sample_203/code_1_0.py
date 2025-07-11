def solve_heesch_number():
  """
  This function determines and prints the Heesch number for the given polygons.

  The polygon in the figure has 180-degree rotational symmetry (point symmetry).
  A theorem in geometry states that any polygon with point symmetry can tile the plane.
  If a polygon can tile the plane, it can be surrounded by an infinite number of
  layers of itself. By definition, its Heesch number is infinity.

  Since all three polygons in the image are identical, the Heesch number for each
  is infinity.
  """

  # The Heesch number for a single polygon that tiles the plane.
  # We use the string "infinity" to represent ∞.
  heesch_number_for_one_polygon = "infinity"

  # There are three identical polygons.
  number_of_polygons = 3

  # Create a list of the results.
  results = [heesch_number_for_one_polygon] * number_of_polygons

  # Print the results separated by commas.
  print(", ".join(results))

solve_heesch_number()