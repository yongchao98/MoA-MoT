import sympy as sp

def main():
    """
    This script symbolically derives the M_32 element of the inertia matrix
    for the given RPR robot.
    """
    # Define symbolic variables
    q1, q2, q3 = sp.symbols('q1 q2 q3')
    L1, dc3, m3 = sp.symbols('L1 d_c3 m3')

    # Define helper variables for trigonometric functions
    c1, s1 = sp.cos(q1), sp.sin(q1)
    c3, s3 = sp.cos(q3), sp.sin(q3)

    # --- Step 1: Define Kinematics ---

    # Transformation from frame 0 to 1
    # R_0_1 based on the frame definitions in the figure
    R_0_1 = sp.Matrix([[-s1, 0, c1],
                       [ c1, 0, s1],
                       [  0, 1,  0]])
    p_0_o1 = sp.Matrix([0, 0, L1])

    # Position of origin of frame 2 in frame 0
    p_1_o2 = sp.Matrix([q2, 0, 0])
    p_0_o2 = p_0_o1 + R_0_1 * p_1_o2

    # Transformation from frame 1 to 2 (constant rotation)
    R_1_2 = sp.Matrix([[0, 0, 1],
                       [1, 0, 0],
                       [0, 1, 0]])

    # Transformation from frame 2 to 3
    R_2_3 = sp.Matrix([[c3, -s3, 0],
                       [s3,  c3, 0],
                       [ 0,   0, 1]])

    # Combined rotation matrix from 0 to 2 and 0 to 3
    R_0_2 = R_0_1 * R_1_2
    R_0_3 = R_0_2 * R_2_3
    
    # --- Step 2: Locate Center of Mass of Link 3 ---
    # Position of CoM of link 3 in frame 3
    p_3_c3 = sp.Matrix([dc3, 0, 0])
    
    # Position of CoM of link 3 in frame 0
    # p_0_c3 = p_0_o2 + R_0_3 * p_3_c3_in_frame_3 is incorrect
    # p_0_c3 = p_0_o3 + R_0_3 * p_3_c3_in_frame_3
    # o3 = o2, so p_0_o3 = p_0_o2
    p_0_c3 = p_0_o2 + R_0_3 * sp.Matrix([1,0,0]) * dc3

    # --- Step 3: Compute Jacobian Columns ---
    # Since joint 2 is prismatic, its motion does not induce angular velocity on subsequent links.
    # Therefore, J_w3_2 (the 2nd column of the angular velocity Jacobian for link 3) is zero.
    # This means the angular part of M_32 is zero. We only need the linear part.

    # J_v_c3_2: Partial derivative of p_0_c3 with respect to q2
    J_v_c3_2 = p_0_c3.diff(q2)
    
    # J_v_c3_3: Partial derivative of p_0_c3 with respect to q3
    J_v_c3_3 = p_0_c3.diff(q3)

    # --- Step 4: Calculate M_32 ---
    # M_32 = m3 * (J_v_c3_3.T * J_v_c3_2) + (angular_term = 0)
    # Note: M is symmetric, so M_32 = M_23. We compute M_23 here.
    M_23_scalar = (J_v_c3_2.T * J_v_c3_3)[0]
    
    # Final expression for M_32
    M_32 = m3 * sp.simplify(M_23_scalar)

    # --- Output the result ---
    print("The expression for the inertia matrix entry M_32 is:")
    print(f"M_32 = {M_32}")
    
    final_result_value = 0
    # The final equation is M_32 = 0
    print("\nThe final equation is:")
    print(f"M_32 = {final_result_value}")


if __name__ == '__main__':
    main()
