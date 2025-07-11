def solve():
  """
  This function solves the problem about the possible number of vertices
  of a convex polyhedron with three quadrilateral projections.
  The possible numbers of vertices for such a polyhedron are 4, 5, 6, and 8.
  The code will print these numbers.
  """
  possible_vertices = [4, 5, 6, 8]
  print("The set of possible numbers of vertices is:")
  for v in possible_vertices:
    print(v)

solve()