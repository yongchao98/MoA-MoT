def solve():
  """
  This function determines and prints the set of possible numbers of vertices
  for a 3D convex polyhedron that can be projected onto three planes in
  general position as a quadrilateral.
  """
  # The set of possible numbers of vertices, V, is composed of:
  # 1. V = 4 (e.g., a tetrahedron)
  # 2. V = 5 (e.g., a square pyramid)
  # 3. All even integers V >= 6 (e.g., n-gonal prisms for n>=3, octahedron for V=6, cube for V=8)
  # It has been shown that odd numbers of vertices V > 5 are not possible.
  
  possible_numbers_description = "The set of possible numbers of vertices is {5} U {V | V is an even integer and V >= 4}."
  
  # For clarity and to list the first few, we can represent it as:
  print("The set of possible numbers of vertices for such a polyhedron P is:")
  print("V = 4")
  print("V = 5")
  print("V = 6")
  print("V = 8")
  print("V = 10")
  print("V = 12")
  print("and so on for all subsequent even numbers.")
  print("\nMore formally, the set is composed of the number 5, and all even integers greater than or equal to 4.")

solve()