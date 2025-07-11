import sympy as sp

def solve_inertia_matrix_entry():
    """
    This function symbolically calculates the M_32 entry of the inertia matrix
    for the RPR robot shown in the figure.
    """
    # Define symbolic variables for joint angles, masses, and dimensions
    q1, q3 = sp.symbols('q1 q3')
    m3 = sp.Symbol('m_3')
    d_c3 = sp.Symbol('d_c3')

    # Define trigonometric functions for convenience
    c1, s1 = sp.cos(q1), sp.sin(q1)
    c3, s3 = sp.cos(q3), sp.sin(q3)

    # Step 1: Define Rotation Matrices based on the provided frames
    # R_0_1: Rotation from frame {0} to {1}
    # x1 axis = [c1, s1, 0], y1 axis = [0, 0, 1], z1 axis = [s1, -c1, 0]
    R_0_1 = sp.Matrix([
        [c1, 0, s1],
        [s1, 0, -c1],
        [0,  1, 0]
    ])

    # R_1_2: Rotation from frame {1} to {2} (constant)
    # x2 axis = y1, y2 axis = x1, z2 axis = -z1
    R_1_2 = sp.Matrix([
        [0, 1,  0],
        [1, 0,  0],
        [0, 0, -1]
    ])

    # R_2_3: Rotation from frame {2} to {3} (around z2 by q3)
    R_2_3 = sp.Matrix([
        [c3, -s3, 0],
        [s3,  c3, 0],
        [0,   0, 1]
    ])

    # Step 2: Calculate the required Jacobian columns for the center of mass of link 3
    # J_vc3_2: Contribution of joint 2 (prismatic along x1) to the linear velocity of CoM 3
    # This is the axis of the prismatic joint in the base frame {0}.
    J_vc3_2 = R_0_1 * sp.Matrix([1, 0, 0])

    # J_vc3_3: Contribution of joint 3 (revolute around z2) to the linear velocity of CoM 3
    # Formula: J_vc3_3 = z2_axis x (p_c3 - p_2)

    # First, find the axis of joint 3 (z2) in the base frame {0}
    R_0_2 = R_0_1 * R_1_2
    z2_axis_0 = R_0_2 * sp.Matrix([0, 0, 1])

    # Second, find the vector from joint 3's origin (p2) to CoM 3 (pc3) in frame {0}
    # The vector in frame {3} is [d_c3, 0, 0]^T
    p_c3_in_frame3 = sp.Matrix([d_c3, 0, 0])
    # We need to rotate this vector to frame {0}
    R_0_3 = R_0_2 * R_2_3
    r_2_c3_0 = R_0_3 * p_c3_in_frame3
    
    # Now, calculate J_vc3_3 using the cross product
    J_vc3_3 = z2_axis_0.cross(r_2_c3_0)

    # Step 3: Calculate M_32 = m_3 * J_vc3_3^T * J_vc3_2
    # The formula for M_ij is symmetric, but using the correct indices is good practice.
    # M_32 = m_3 * J_vc3,3.T * J_vc3,2
    M_32_expr = m3 * (J_vc3_3.transpose() * J_vc3_2)

    # The result is a 1x1 matrix, so we extract the element
    M_32_final = sp.simplify(M_32_expr[0])

    # Step 4: Print the final expression and its components
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    # Using sp.pretty_print for better mathematical formatting
    sp.pprint(M_32_final, use_unicode=True)
    
    print("\n--- Equation Breakdown ---")
    print(f"Final Equation: M_32 = {M_32_final}")
    print("This equation shows the dynamic coupling between the prismatic joint 2 and the revolute joint 3.")
    print("The components are:")
    print(f"- {m3}: The mass of link 3.")
    print(f"- {d_c3}: The distance from the axis of joint 3 to the center of mass of link 3.")
    print(f"- cos({q3}): The cosine of the angle of joint 3, which varies with the robot's configuration.")
    print("The numerical coefficient of this expression is 1.")

if __name__ == '__main__':
    solve_inertia_matrix_entry()
