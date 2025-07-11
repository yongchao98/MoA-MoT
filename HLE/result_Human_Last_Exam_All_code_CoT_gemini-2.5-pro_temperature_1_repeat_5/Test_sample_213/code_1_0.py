# The user wants the original, unmodified line of code.
# Based on the analysis of standard 3D rotation matrices, the 'Rz' matrix is missing a negative sign.
# The buggy line is:
# mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(sin(a.z),cos(a.z),0),vec3(0,0,1));
# The second vector should be vec3(-sin(a.z), cos(a.z), 0) for a correct Z-axis rotation.
# Therefore, the original line had a '-' character before 'sin(a.z)'.

original_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));"
print(original_line)