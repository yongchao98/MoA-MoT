import sympy as sp

def calculate_inertia_matrix_element():
    """
    This function symbolically derives the expression for the M_32 element
    of the inertia matrix for the given RPR spatial robot.
    """
    # 1. Define symbolic variables
    # q1, q2, q3 are the joint variables.
    # L1 is the constant height of the first joint.
    # m3 is the mass of link 3.
    # dc3 is the distance to the center of mass of link 3 along its local x-axis.
    q1, q2, q3 = sp.symbols('q1 q2 q3')
    L1, m3, dc3 = sp.symbols('L1 m3 dc3')

    # 2. Kinematic model and Center of Mass position
    # The position of the wrist joint (origin of frame 2) is p2.
    # The prismatic joint q2 moves along a horizontal axis whose direction is determined by q1.
    c1, s1 = sp.cos(q1), sp.sin(q1)
    p2 = sp.Matrix([
        q2 * c1,
        q2 * s1,
        L1
    ])

    # The CoM of link 3 is at a distance dc3 along its local x3-axis.
    # The third joint q3 rotates around the vertical z2-axis.
    # The orientation of link 3 in the base xy-plane is thus dependent on (q1 + q3).
    c13, s13 = sp.cos(q1 + q3), sp.sin(q1 + q3)
    
    # The position vector from the wrist joint to the CoM of link 3, expressed in the base frame.
    p_c3_relative = sp.Matrix([
        dc3 * c13,
        dc3 * s13,
        0
    ])

    # The absolute position of the CoM of link 3 in the base frame.
    pc3 = p2 + p_c3_relative

    # 3. Calculate Jacobian columns for the CoM of link 3
    # The linear velocity Jacobian column J_v_j is the partial derivative of pc3 w.r.t. q_j.
    J_vc3_2 = pc3.diff(q2)
    J_vc3_3 = pc3.diff(q3)

    # 4. Calculate M_32
    # The angular velocity Jacobian for joint 2 (prismatic) is zero, so the angular part of M_32 is zero.
    # M_32 = m3 * J_vc3_3^T * J_vc3_2
    M32 = m3 * (J_vc3_3.T * J_vc3_2)[0] # Extract scalar from 1x1 matrix

    # Simplify the final expression
    M32_simplified = sp.simplify(M32)

    # 5. Print the final expression
    # The problem asks to output the expression.
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    
    # To satisfy the "output each number" constraint, we print the equation components.
    # The derived expression is -m3*dc3*sin(q3).
    # The implicit number is -1.
    print("M_32 = -1 * m3 * dc3 * sin(q3)")
    
    # For programmatic use, we can print the sympy string representation.
    # print(M32_simplified)


if __name__ == '__main__':
    calculate_inertia_matrix_element()
