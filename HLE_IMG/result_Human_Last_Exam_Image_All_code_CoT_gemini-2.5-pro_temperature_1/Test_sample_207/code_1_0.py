import sympy

def get_m32_expression():
    """
    Calculates and prints the symbolic expression for the M_32 entry of the
    inertia matrix for the given RPR robot.

    The derivation is based on the formula:
    M_jk = sum over i [ m_i * J_v_ci,j^T * J_v_ci,k + J_omega_i,j^T * R_i * I_ci * R_i^T * J_omega_i,k ]

    For M_32, the sum is only for link i=3.
    The rotational term is zero because joint 2 is prismatic, so its corresponding
    angular velocity Jacobian column (J_omega_3,2) is zero.
    The translational term is derived from the kinematics based on the provided figure.
    The final result depends only on the properties of link 3 and the joint angle q_3.
    """
    # Define symbols for the variables in the expression
    m_3 = sympy.Symbol('m_3')
    d_c3 = sympy.Symbol('d_c3')
    q_3 = sympy.Symbol('q_3')

    # The derived expression for M_32
    M_32 = -m_3 * d_c3 * sympy.sin(q_3)

    # Print the final expression
    print("The expression for the entry M_32 of the robot inertia matrix is:")
    print(f"M_32 = {M_32}")

    # As requested, output each component of the final equation
    print("\nThe final equation can be written as:")
    print(f"-1 * {m_3} * {d_c3} * sin({q_3})")

# Execute the function to print the result
get_m32_expression()