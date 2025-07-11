# The user wants me to identify and correct a single-character bug in a GLSL shader.
# The bug is in the definition of the rotation matrix for the Z-axis, Rz.
# The original code's Rz is not an orthogonal rotation matrix because it's missing a negative sign.
# By restoring the negative sign, the matrix becomes a standard rotation matrix.

# The buggy line is:
# mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(sin(a.z),cos(a.z),0),vec3(0,0,1));
# The corrected line, with the single missing character restored, is:
# mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));

# The user asked for the original, unmodified line to be printed exactly.
print('mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));')