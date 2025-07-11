# The user wants the original, unmodified line of code.
# Based on the analysis, the bug is a missing '-' sign in the definition of the Rz rotation matrix.
# The matrix definition is split across two lines in the original prompt, but it's a single statement.
# The prompt asks for the original, unmodified line exactly.
# The buggy line is: mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(sin(a.z),cos(a.z),0),vec3(0,0,1));
# The corrected line adds a '-' before the second sin(a.z).

original_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));"
print(original_line)