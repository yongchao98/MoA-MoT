def solve():
  """
  This function prints the corrected line of GLSL code.
  The bug was a missing '-' sign in the Z-axis rotation matrix definition,
  which is crucial for a correct rotation transformation.
  """
  correct_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));"
  print(correct_line)

solve()