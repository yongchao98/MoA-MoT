def solve():
    """
    This function prints the corrected line of GLSL code.
    The bug is in the definition of the Z-axis rotation matrix 'Rz'.
    In GLSL, mat3 constructors take column vectors. The standard column-major
    rotation matrix for the Z-axis is:
    [cos(a)  -sin(a)  0]
    [sin(a)   cos(a)  0]
    [0        0       1]

    This corresponds to the arguments:
    vec3(cos(a), sin(a), 0),
    vec3(-sin(a), cos(a), 0),
    vec3(0, 0, 1)

    The provided code is missing the negative sign '-' in the second vector,
    which constitutes the single-character bug.
    """
    original_line = "mat3 Rz=mat3(vec3(cos(a.z),sin(a.z),0),vec3(-sin(a.z),cos(a.z),0),vec3(0,0,1));"
    print(original_line)

solve()