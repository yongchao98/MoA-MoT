# The corrected line fixes the Z-axis rotation matrix by adding the missing negative sign.
# Without this sign, the matrix does not perform a pure rotation, leading to visual bugs.
correct_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));"
print(correct_line)