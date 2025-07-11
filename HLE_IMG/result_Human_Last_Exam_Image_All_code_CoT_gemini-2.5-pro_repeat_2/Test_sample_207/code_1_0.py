import sympy
from sympy import sin, cos, symbols, Matrix

def solve_inertia_matrix_element():
    """
    This function symbolically derives the M_32 element of the inertia matrix
    for the RPR robot shown in the figure.
    """
    # 1. Define symbolic variables based on the problem description.
    # m_3: mass of link 3
    # d_c3: distance of the Center of Mass (CoM) of link 3 from the origin of frame 3, along the x3-axis.
    # q_1, q_3: revolute joint angles
    # q_2, d_1: prismatic joint variable and vertical offset (these will cancel out in the final M_32 expression)
    m_3, d_c3 = symbols('m_3 d_c3')
    q_1, q_2, q_3 = symbols('q_1 q_2 q_3')
    d_1 = symbols('d_1') 

    c1, s1 = cos(q_1), sin(q_1)
    c3, s3 = cos(q_3), sin(q_3)

    # 2. Define key kinematic vectors and matrices in the base frame {0}.
    # k2: Direction vector of prismatic joint 2 (along x1 axis).
    k2 = Matrix([c1, s1, 0])

    # k3: Axis vector of revolute joint 3 (along z2 axis).
    k3 = Matrix([s1, -c1, 0])
    
    # p2: Position of the origin of joint 3's frame.
    p2 = Matrix([q_2 * c1, q_2 * s1, d_1])

    # To find the position of the CoM of link 3 (pc3), we need the rotation matrix R_30.
    # The frames are defined such that y1 || z0, and x1, z1 rotate with q1.
    R_10 = Matrix([[c1, 0, s1], [s1, 0, -c1], [0, 1, 0]])
    # Frame 2 is translated from frame 1 but not rotated (R_21 = Identity).
    # Frame 3 is rotated about the z2 axis by q3.
    R_32 = Matrix([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
    R_30 = R_10 * R_32
    
    # Position of CoM of link 3 in its own frame {3}.
    rc3_in_frame3 = Matrix([d_c3, 0, 0])
    
    # Position of CoM of link 3 in the base frame {0}.
    pc3 = p2 + R_30 * rc3_in_frame3
    
    # 3. Calculate the required Jacobian columns for the CoM of link 3.
    # For joint 2 (prismatic):
    J_v32 = k2
    J_w32 = Matrix([0, 0, 0])

    # For joint 3 (revolute):
    # Vector from joint 3 axis to the CoM of link 3.
    r_from_j3_to_c3 = pc3 - p2
    J_v33 = k3.cross(r_from_j3_to_c3)
    J_w33 = k3

    # 4. Calculate M_32 using the general formula for inertia matrix elements.
    # M_32 = m_3 * J_v32^T * J_v33 + J_w32^T * R_3 * I_c3 * R_3^T * J_w33
    # The angular term is zero because J_w32 is the zero vector.
    M32_linear_term = m_3 * (J_v32.T * J_v33)[0] # Extract the scalar from the 1x1 matrix

    # 5. Simplify the expression and print the result.
    M_32 = sympy.simplify(M32_linear_term)

    # Print the final expression in a formatted way.
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    # The output format follows the structure of the equation.
    print(f"M_32 = {M_32}")
    
solve_inertia_matrix_element()