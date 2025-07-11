def solve():
    """
    This function prints the original, unmodified line of GLSL code.
    The bug was a single missing '-' character in the Z-axis rotation matrix definition.
    """
    original_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));"
    print(original_line)

solve()