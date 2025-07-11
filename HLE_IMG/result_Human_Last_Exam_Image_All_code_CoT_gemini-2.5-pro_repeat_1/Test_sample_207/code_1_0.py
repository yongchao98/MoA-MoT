import sympy as sp

def solve_inertia_matrix_element():
    """
    This function calculates and prints the expression for the M_32 element
    of the inertia matrix for the given RPR robot.

    The derivation is based on Jacobian formulation of robot dynamics. The M_jk
    element of the inertia matrix M(q) is given by:
    M_jk = Sum_i [ m_i * J_v_ci(j)^T * J_v_ci(k) + J_w_i(j)^T * R_i * I_ci * R_i^T * J_w_i(k) ]

    For M_32, j=3 and k=2. Since joint 2 is prismatic, its corresponding
    angular velocity Jacobian column, J_w(2), is zero. This simplifies the
    expression significantly, as the second term in the summation becomes zero.

    Furthermore, only links affected by both joints 2 and 3 will contribute to
    the sum. In this RPR manipulator, only link 3 is affected by both.
    Therefore, the expression simplifies to:
    M_32 = m_3 * J_v_c3(3)^T * J_v_c3(2)

    After deriving the kinematic transformations from the figure and calculating
    the required Jacobians, the final result is obtained.
    """

    # Define the symbolic variables used in the expression
    # m_3: mass of link 3
    # d_c3: distance of the center of mass of link 3 from the origin of frame 3, along the x3 axis
    # q_3: joint angle for the third (revolute) joint
    m3 = sp.Symbol('m_3')
    dc3 = sp.Symbol('d_c3')
    q3 = sp.Symbol('q_3')

    # The derived expression for M_32
    # The calculation involves finding the Jacobians for the center of mass of link 3
    # and taking their dot product. The result simplifies to the following expression.
    M_32 = m3 * dc3 * sp.cos(q3)

    # Print the final result in a formatted equation
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    
    # We construct and print the equation M_32 = expression
    # The numbers in the equation are implicitly 1.
    final_equation = sp.Eq(sp.Symbol('M_32'), M_32, evaluate=False)
    
    # Pretty print for better readability
    sp.pprint(final_equation, use_unicode=True)
    
    # For clarity, also print the formula in a simple text format,
    # explicitly showing all components of the final equation as requested.
    print("\nIn standard notation, the equation is:")
    print(f"M_32 = 1 * {m3} * {dc3} * cos({q3})")


if __name__ == "__main__":
    solve_inertia_matrix_element()