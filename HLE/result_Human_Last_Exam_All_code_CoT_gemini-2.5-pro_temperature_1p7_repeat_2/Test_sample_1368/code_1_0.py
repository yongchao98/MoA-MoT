def solve():
  """
  This function prints the shortest possible sequence of commands to construct a square
  inscribed in a circle.
  """
  # L: Draw a Line through the center O and point A on the circumference to find the opposite vertex C.
  # C: Draw a Circle centered at A with radius AC.
  # C: Draw a Circle centered at C with radius CA.
  # L: Draw a Line through the two intersection points of the previous two circles.
  # This line is the second diagonal, giving the remaining two vertices.
  construction_sequence = "LCCL"
  print(construction_sequence)

solve()