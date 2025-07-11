import sympy
from sympy import symbols, cos, sin, Matrix, simplify

def solve_m32_expression():
    """
    This function symbolically derives and prints the expression for the M_32
    element of the inertia matrix for the given RPR robot.
    """
    # Define symbolic variables for masses, lengths, and joint variables
    # m3: mass of link 3
    # l3: length of link 3 (DH parameter a_2)
    # dc3: distance to the center of mass of link 3 along its x-axis (x_3)
    # q1, q3: joint angles for joint 1 and 3
    m3, l3, dc3 = symbols('m_3 l_3 d_c3')
    q1, q3 = symbols('q_1 q_3')

    # Define sines and cosines for readability
    s1, c1 = sin(q1), cos(q1)
    s3, c3 = sin(q3), cos(q3)

    # --- Derivation based on the Jacobian method ---

    # M_32 is determined solely by the properties of link 3.
    # The angular velocity contribution is zero because joint 2 is prismatic.
    # M_32 = m_3 * J_v_c3,2^T * J_v_c3,3

    # J_v_c3,2 is the axis of prismatic joint 2, which is z_1.
    # From the DH parameter analysis of the robot:
    # z_1 is the z-axis of frame 0 rotated by q1 about z0, then by pi/2 about x1.
    # This results in z_1 = [sin(q1), -cos(q1), 0]^T
    z_1 = Matrix([s1, -c1, 0])
    J_v_c3_2 = z_1

    # J_v_c3,3 = z_2 x (p_c3 - p_2), where z_2 is the axis of revolute joint 3.
    # z_2 = [cos(q1), sin(q1), 0]^T
    z_2 = Matrix([c1, s1, 0])

    # p_c3 - p_2 is the vector from the origin of frame 2 to the CoM of link 3,
    # expressed in the base frame.
    # First, find this vector in frame 2 coordinates (r_c3_in_2):
    # The CoM is at [dc3, 0, 0] in frame 3. The origin of frame 3 is at [l3, 0, 0] in frame 2.
    # r_c3_in_2 = [l3 + dc3*cos(q3), dc3*sin(q3), 0]^T
    # Then, rotate this vector to the base frame using R_0^2:
    R_0_2 = Matrix([[0, s1, c1], [0, -c1, s1], [1, 0, 0]])
    r_c3_in_2 = Matrix([l3 + dc3*c3, dc3*s3, 0])
    p_c3_minus_p_2 = R_0_2 * r_c3_in_2

    # Now, calculate the cross product for J_v_c3,3
    J_v_c3_3 = z_2.cross(p_c3_minus_p_2)

    # Finally, compute M_32 by the dot product
    M_32 = m3 * J_v_c3_2.dot(J_v_c3_3)

    # Simplify the final expression
    M_32_simplified = simplify(M_32)

    # Print the final expression for M_32
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    sympy.pprint(M_32_simplified)
    # To explicitly show the final equation form:
    print("\nIn equation form:")
    print(f"M_32 = {M_32_simplified}")


solve_m32_expression()