import sympy as sp

def provide_inertia_matrix_expression():
    """
    This function provides the symbolic expression for the M_32 entry of the
    inertia matrix for the given RPR robot, based on a detailed kinematic and
    dynamic analysis.

    Derivation Summary:
    1. The inertia matrix entry M_ij is given by the sum over links k (from max(i, j) to n):
       M_ij = m_k * J_v_ki^T * J_v_kj + J_w_ki^T * R_k * I_k * R_k^T * J_w_kj

    2. For M_32, we only need to consider link 3 (k=3). Since joint 2 is prismatic, its
       contribution to the angular velocity of link 3 is zero. Thus, the angular velocity
       Jacobian J_w_3,2 is the zero vector, and the second term of the formula vanishes.
       The expression simplifies to: M_32 = m_3 * J_v_c3,2^T * J_v_c3,3

    3. Jacobian for joint 2 (prismatic): The motion is along axis y1, which is parallel to
       the world z0 axis. So, J_v_c3,2 = [0, 0, 1]^T.

    4. Jacobian for joint 3 (revolute): This is found using the formula J_v_c3,3 = z_2 x (p_c3 - p_2).
       By defining the coordinate frames as shown in the figure and performing the vector
       cross product, this Jacobian is calculated.

    5. Dot Product: The final expression is obtained by computing the dot product of the
       two Jacobians and multiplying by the mass of link 3, m_3.

    The resulting expression is M_32 = m_3 * d_c3 * cos(q_3).
    """

    # Define symbols for the parameters in the final equation.
    # m_3: mass of link 3
    # d_c3: distance to the center of mass of link 3 along its x3 axis
    # q_3: joint variable for the third (revolute) joint
    m_3 = sp.Symbol("m_3")
    d_c3 = sp.Symbol("d_c3")
    q_3 = sp.Symbol("q_3")

    # The final expression for M_32, as derived from the analysis.
    expression = m_3 * d_c3 * sp.cos(q_3)

    # Print the final result in a clear, readable format, showing each component.
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    print(f"M_32 = {m_3} * {d_c3} * cos({q_3})")


# Run the function to display the answer.
provide_inertia_matrix_expression()