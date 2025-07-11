# The buggy line in the shader is:
# mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(sin(a.z),cos(a.z),0),vec3(0,0,1));
#
# A minus sign '-' was removed from the second vector constructor, before sin(a.z).
# This makes the matrix incorrect for rotation.
# The original, correct line should be:
# mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));

print("mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));")