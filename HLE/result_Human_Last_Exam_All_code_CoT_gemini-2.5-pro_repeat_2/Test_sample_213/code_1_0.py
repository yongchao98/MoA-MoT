def solve():
  """
  This function identifies and corrects the buggy line in the GLSL code.
  The bug is a missing '-' sign in the Z-axis rotation matrix definition.
  """
  
  # The buggy line provided by the user.
  # buggy_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(sin(a.z),cos(a.z),0),vec3(0,0,1));"
  
  # The original, correct line with the '-' character restored.
  correct_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));"
  
  print(correct_line)

solve()